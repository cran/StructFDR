\name{StructFDR}
\alias{StructFDR}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
False Discovery Rate (FDR) Control Integrating a General Prior Structure
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The procedure is based on an empirical Bayes hierarchical model, where a structure-based prior distribution is designed to utilize the prior structure information. A moderated
statistic based on posterior mean is used for permutation-based FDR control. By borrowing information from neighboring features defined based a distance metric, it is able to 
improve the statistical power of detecting associated features while controlling the FDR at desired levels. 
}
\usage{
StructFDR(X, Y, D, test.func, perm.func, eff.sign = TRUE, B = 20, c.cutoff = 0.001, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
a data matrix, rows are the features and columns are the samples.
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
a vector of the phenotypic values, where association tests are being assessed
}
  \item{D}{
%%     ~~Describe \code{tree} here~~
a matrix of pairwise distances between features. It determines how the features are related and the neighborhood for information borrowing.
}
  \item{test.func}{
%%     ~~Describe \code{test.func} here~~
a function that performs the actual tests. It takes \code{X}, \code{Y} and ... as the inputs, and returns a list with two slots \code{p.value} and \code{e.sign}, which are vectors of p-values 
and signs of the effects.
}
  \item{perm.func}{
%%     ~~Describe \code{perm.func} here~~
a function that performs the permutation tests. It takes \code{X}, \code{Y} and ... as the inputs, and returns a list with two slots \code{X} and \code{Y}, which contain the permuted data.
}
  \item{eff.sign}{
%%     ~~Describe \code{eff.sign} here~~
a logical value indicating whether the direction of the effects should be considered. If it is true (default), negative and positive effects provide conflicting 
information. 
}
  \item{B}{
%%     ~~Describe \code{B} here~~
number of the permutations.  The default is 20.
}
  \item{c.cutoff}{
%%     ~~Describe \code{c.cutoff} here~~
a threshold to which the standard deviation of the moderated statistic is compared. If the standard deviation is smaller than the threshold,
it indicates potential "over-adjustment". The standard BH procedure will be used.
}

  \item{...}{
%%     ~~Describe \code{c.cutoff} here~~
further arguments such as covariates to be passed to \code{test.func}
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{p.adj }{StructFDR adjusted p-values.}
\item{p.unadj }{raw p-values.}
\item{z.adj }{moderated z-values. The scale may be different from the raw z-values.}
\item{z.unadj }{raw z-values.}
\item{k, rho }{the estimates of the hyperparameters.  The magnitude indicates the informativeness of the prior structure.}

}

\references{
%% ~put references to the literature/web site here ~
Jian Xiao, Hongyuan Cao and Jun Chen (2016). False discovery rate control incorporating phylogenetic tree increases detection power in microbiome-wide multiple testing. Submitted.
}
\author{
%%  ~~who you are~~
Jun Chen
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{TreeFDR}}
}
\examples{
require(ape) 
require(nlme)
require(cluster)
require(StructFDR)

# Generate a caelescence tree and partition into 10 clusters
set.seed(1234)
n <- 20
p <- 200
tree <- rcoal(p)
# Pairwise distance matrix. 
D <- cophenetic(tree)
clustering <- pam(D, k=10)$clustering

# Simulate case-control data, assuming cluster 2 is differential  
X.control <- matrix(rnorm(n*p), p, n)
X.case <- matrix(rnorm(n*p), p, n)
eff.size <- rnorm(sum(clustering == 2), 0.5, 0.2)     
X.case[clustering == 2, ] <- X.case[clustering == 2, ] + eff.size
X <- cbind(X.control, X.case)
Y <- gl(2, n) 

# Define testing and permutation function
test.func <- function (X, Y) {
	obj <- apply(X, 1, function(x) {
				ttest.obj <- t.test(x ~ Y)
				c(ttest.obj$p.value, sign(ttest.obj$statistic))
			})
    return(list(p.value=obj[1, ], e.sign=obj[2, ]))
}

perm.func <- function (X, Y) {
	return(list(X=X, Y=sample(Y)))
}

# Call StructFDR
tree.fdr.obj <- StructFDR(X, Y, D, test.func, perm.func)

# Compare StructFDR and BH
tree.fdr.obj$p.adj
tree.fdr.obj$p.adj[clustering == 2]
BH.p.adj <- p.adjust(tree.fdr.obj$p.unadj, 'fdr')
BH.p.adj[clustering == 2]

# Adjusted statistics vs clustering
par(mfrow=c(1, 2))
plot(clustering, tree.fdr.obj$z.unadj)
plot(clustering, tree.fdr.obj$z.adj)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Multiple testing }
\keyword{ False discovery rate}
\keyword{ Genetics } 
\keyword{ Metagenomics }
\keyword{ Genomics }% __ONLY ONE__ keyword per line

