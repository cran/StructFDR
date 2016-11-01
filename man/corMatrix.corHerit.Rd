\name{corMatrix.corHerit}
\alias{corMatrix.corHerit}
\title{Matrix of a corHerit Object
}
\description{This function provides actual computation of the correlation matrix and also extracts the correlation matrix. }
\usage{
\method{corMatrix}{corHerit}(object, covariate = getCovariate(object), ...)
}
\arguments{
  \item{object}{
an object inheriting from class \code{corHerit} representing a correlation structure.
}
  \item{covariate}{
an optional covariate vector (matrix), or list of covariate vectors (matrices), at which values the correlation matrix, or list of correlation matrices, are to be evaluated. Defaults to getCovariate(object).
}
  \item{\dots}{
some methods for this generic require additional arguments. 
}
}
\value{
If covariate is a vector (matrix), the returned value will be an array with the corresponding correlation matrix. If the covariate is a list of vectors (matrices), the returned value will be a list with the correlation matrices corresponding to each component of covariate.
}
\references{Pinheiro, J.C., and Bates, D.M. (2000) "Mixed-Effects Models in S and S-PLUS", Springer.
}


