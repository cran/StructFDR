require(nlme)

# Transform to Z score 
Ztransform <- function(p.value, e.sign, eff.sign=TRUE, tol=1E-15) {
    p.value[p.value <= tol] <- tol
	p.value[p.value >= 1 - tol] <- 1 - tol
	if (eff.sign == TRUE) {
		e.sign[e.sign == 0] <- sample(c(-1, 1), sum(e.sign == 0), replace=T)
		z1 <- qnorm(p.value / 2)
		z2 <- qnorm(1 - p.value / 2)
		z <- ifelse(e.sign > 0, z2, z1)
	} else {
		z <- qnorm(1 - p.value)
	}
	return(z)
}

# Internal functions
# For gls with generic correlation structure
corHerit <- function(value, paras, form = ~1, fixed = TRUE) {
	# Place holder - check the validity of the parameter
	object <- value
	attr(object, "formula") <- form
	attr(object, "fixed") <- fixed
	attr(object, "paras") <- paras
	class(object) <- c("corHerit", "corStruct")
	object
}

Initialize.corHerit <- function (object, data, ...) {	
	# Place holder - check the validity of the parameter
	form <- formula(object)
	if (!is.null(getGroupsFormula(form))) {
		attr(object, "groups") <- getGroups(object, form, data = data)
		attr(object, "Dim") <- Dim(object, attr(object, "groups"))
	} else {
		attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))))
	}
	attr(object, "covariate") <- getCovariate(object, data = data)
	object
}


corMatrix.corHerit <- function (object, covariate = getCovariate(object), ...) {
	
	paras <- attr(object, "paras")
	p <- paras[['p']]
	I <- diag(p)
	hyper <- as.vector(object)
	
	lambda <- 1 / (1 + exp(hyper[1]))
	V <- exp(-((paras[['D']])^2) * (exp(hyper[2])))
	cor.m <- (1 - lambda) * V +  lambda * I
	cor.m
}

# Extract the coefficient
coef.corHerit <- function (object,  ...) {
	
	paras <- attr(object, "paras")
	coefs <- as.vector(object)	
	coef1 <- 1 / (1 + exp(coefs[1]))
	coef1 <- coef1 / (1 - coef1)
	
	coef2 <-  exp(coefs[2])
	coefs <- c(coef1, coef2)	
	names(coefs) <- paste("Hyper", 1:length(coefs), sep="")
	coefs
}


# Estimating the hyper parameters
EstHyper <- function (y, D, init.val=c(0, 0)) {


	obj.gls <- gls(model = y ~ 1, 
			correlation = corHerit(value=init.val, paras = list(p=length(y), D = D)))
	cc <- c(coef(obj.gls$modelStruct$corStruct), obj.gls$coefficients, obj.gls$logLik)
	cc	
}	


AdjStats <- function (y, V, k, mu, fudge=0.005) {
	p <- nrow(V)
	# Add small fudge
	V.inv <- solve(V + fudge * diag(p))
	I <- diag(p)
	y.adj  <- solve(I + k * V.inv) %*% (k * mu * rowSums(V.inv) + y)
	y.adj[, 1]
}

# For microbiome-applications
TreeFDR <- function(X, Y, tree, test.func, perm.func, eff.sign=TRUE, B=20, c.cutoff=1e-3, ...) {
	
	# Patristic distance
    D <- cophenetic(tree)
	
	if (nrow(D) != nrow(X)) {
		stop('The tree is not consistent with the data! Please check! \n')
	}
	
	cat('Perform testing ... \n')
	
	test.obs <- test.func(X, Y, ...)
	
	if (!is.list(test.obs) | !all(c('e.sign', 'p.value') %in% names(test.obs))) {
		stop("test.func should return a list with names e.sign and p.value! Please check!\n")
	}
	
	z.obs <- Ztransform(test.obs$p.value, test.obs$e.sign, eff.sign)

	z.perm <- list()
	for (i in 1:B) {
		perm.obj <- perm.func(X, Y, ...)
		
		if (!is.list(perm.obj) | !all(c('X', 'Y') %in% names(perm.obj))) {
			stop("perm.func should return a list with names X and Y! Please check!\n")
		}
		X.perm <- perm.obj$X
		Y.perm <- perm.obj$Y
		test.perm <- test.func(X.perm, Y.perm, ...)
		z.perm[[i]] <- Ztransform(test.perm$p.value, test.perm$e.sign, eff.sign)
	}
	
	# Estimate the hyper parameters.
	cat('Estimating hyperparameter ... \n')
	error <- try(
			obj <- EstHyper(y=z.obs, D=D)
	)	
	
	# If error happens during estimating parameters, TreeFDR method reduces to BH method.
	if (inherits(error, "try-error")) {
		warning('Hyperparameter estimation failed! B-H procedure will be used!\n')
		p.adj <- p.adjust(test.obs$p.value,'fdr')
		k <- NULL
		rho <- NULL
	} else {
		cat('Adjustment and Permutation test ...\n')
		k <- obj[1]
		rho <- obj[2]
		mu <- obj[3]	
		V <- exp(-1 * rho *D^2)
		
		stat.o <- AdjStats(y=z.obs, V=V, k=k, mu=mu)
		z.adj <- stat.o
		if (sd(stat.o) <= c.cutoff) {
			warning('Possible over-adjustment! B-H procedure will be used!\n')
			p.adj <- p.adjust(test.obs$p.value,'fdr')
		} else {
			stat.p <- lapply(1:B, function (i) AdjStats(y=z.perm[[i]], V=V, k=k, mu=mu))	
			
			if (eff.sign == TRUE) {
				stat.o <- abs(stat.o)
				stat.p <- lapply(stat.p, abs)
			}
			
			ord <- order(stat.o, decreasing=T)
			stat.o <- stat.o[ord]
			p.adj <-  sapply(stat.o, function(x) sum(unlist(stat.p) >= x)/length(stat.p)) / (1:length(stat.o)) 
			
			# Conservative
			# p.adj <- pmin(1, cummax(p.adj))[order(ord)]
			p.adj <- pmin(1, rev(cummin(rev(p.adj))))[order(ord)]
		}
		
	}
	cat('Done!\n')
	return(list(p.adj=p.adj, p.unadj=test.obs$p.value, z.adj=z.adj, z.unadj=z.obs, k=k, rho=rho))	
	
}


# For other genomics applications
StructFDR <- function(X, Y, D, test.func, perm.func, eff.sign=TRUE, B=20, c.cutoff=1e-3, ...) {
	
	if (nrow(D) != nrow(X)) {
		stop("The distance matrix is not consistent with the data! Please check! \n")
	}
	
	cat('Perform testing ... \n')
	test.obs <- test.func(X, Y, ...)

	if (!is.list(test.obs) | !all(c('e.sign', 'p.value') %in% names(test.obs))) {
		stop("test.func should return a list with names e.sign and p.value! Please check!\n")
	}
	
	z.obs <- Ztransform(test.obs$p.value, test.obs$e.sign, eff.sign)
	
	z.perm <- list()
	for (i in 1:B) {
		perm.obj <- perm.func(X, Y, ...)
		
		if (!is.list(perm.obj) | !all(c('X', 'Y') %in% names(perm.obj))) {
			stop("perm.func should return a list with names X and Y! Please check!\n")
		}
		
		X.perm <- perm.obj$X
		Y.perm <- perm.obj$Y
		test.perm <- test.func(X.perm, Y.perm, ...)
		z.perm[[i]] <- Ztransform(test.perm$p.value, test.perm$e.sign, eff.sign)
	}
	
	# Estimate the hyper parameters.
	cat('Estimating hyperparameter ... \n')
	error <- try(
			obj <- EstHyper(y=z.obs, D=D)
	)	
	
	# If error happens during estimating parameters, TreeFDR method reduces to BH method.
	if (inherits(error, "try-error")) {
		warning('Hyperparameter estimation failed! B-H procedure will be used!\n')
		p.adj <- p.adjust(test.obs$p.value,'fdr')
		k <- NULL
		rho <- NULL
	} else {
		cat('Adjustment and Permutation test ...\n')
		k <- obj[1]
		rho <- obj[2]
		mu <- obj[3]	
		
		V <- exp(-1 * rho *D^2)
		
		stat.o <- AdjStats(y=z.obs, V=V, k=k, mu=mu)
		z.adj <- stat.o
		if (sd(stat.o) <= c.cutoff) {
			warning('Possible over-adjustment! B-H procedure will be used!\n')
			p.adj <- p.adjust(test.obs$p.value,'fdr')
		} else {
			stat.p <- lapply(1:B, function (i) AdjStats(y=z.perm[[i]], V=V, k=k, mu=mu))	
			
			if (eff.sign == TRUE) {
				stat.o <- abs(stat.o)
				stat.p <- lapply(stat.p, abs)
			} 
			
			ord <- order(stat.o, decreasing=T)
			stat.o <- stat.o[ord]
			p.adj <-  sapply(stat.o, function(x) sum(unlist(stat.p) >= x)/length(stat.p)) / (1:length(stat.o)) 
			
			# Conservative
			# p.adj <- pmin(1, cummax(p.adj))[order(ord)]
			p.adj <- pmin(1, rev(cummin(rev(p.adj))))[order(ord)]
		}
		
	}
	cat('Done!\n')
	return(list(p.adj=p.adj, p.unadj=test.obs$p.value, z.adj=z.adj, z.unadj=z.obs, k=k, rho=rho))	
}


## Example
#require(ape) 
#require(nlme)
#require(cluster)
#
## Generate a caelescence tree and partition into 10 clusters
#set.seed(1234)
#n <- 20
#p <- 200
#tree <- rcoal(p)
#D <- cophenetic(tree)
#clustering <- pam(D, k=10)$clustering
#
## Simulate case-control data, assuming cluster 2 is differential  
#X.control <- matrix(rnorm(n*p), p, n)
#X.case <- matrix(rnorm(n*p), p, n)
#eff.size <- rnorm(sum(clustering == 2), 0.5, 0.2)     
#X.case[clustering == 2, ] <- X.case[clustering == 2, ] + eff.size
#X <- cbind(X.control, X.case)
#Y <- gl(2, n) 
#
## Define testing and permutation function function
#test.func <- function (X, Y) {
#	obj <- apply(X, 1, function(x) {
#				ttest.obj <- t.test(x ~ Y)
#				c(ttest.obj$p.value, sign(ttest.obj$statistic))
#			})
#    return(list(p.value=obj[1, ], e.sign=obj[2, ]))
#}
#
#perm.func <- function (X, Y) {
#	return(list(X=X, Y=sample(Y)))
#}
#
## Call TreeFDR
#tree.fdr.obj <- TreeFDR(X, Y, tree, test.func, perm.func)
#
## Compare TreeFDR and BH
#tree.fdr.obj$p.adj
#tree.fdr.obj$p.adj[clustering == 2]
#BH.p.adj <- p.adjust(tree.fdr.obj$p.unadj, 'fdr')
#BH.p.adj[clustering == 2]
#
## Adjusted statistics vs clustering
#par(mfrow=c(1, 2))
#plot(clustering, tree.fdr.obj$z.unadj)
#plot(clustering, tree.fdr.obj$z.adj)
#
## Call StructFDR
#tree.fdr.obj <- StructFDR(X, Y, D, test.func, perm.func)
#
## Compare TreeFDR and BH
#tree.fdr.obj$p.adj
#tree.fdr.obj$p.adj[clustering == 2]
#BH.p.adj <- p.adjust(tree.fdr.obj$p.unadj, 'fdr')
#BH.p.adj[clustering == 2]
#
## Adjusted statistics vs clustering
#par(mfrow=c(1, 2))
#plot(clustering, tree.fdr.obj$z.unadj)
#plot(clustering, tree.fdr.obj$z.adj)
