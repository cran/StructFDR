\name{Initialize.corHerit}
\alias{Initialize.corHerit}
\title{Initialize corHerit Object
}
\description{
This method initializes object by evaluating its associated covariate(s) and grouping factor, if any is present, in data, calculating various dimensions and constants used by optimization algorithms involving corStruct objects (see the appropriate Dim method documentation), and assigning initial values for the coefficients in object, if none were present.}
\usage{
## S3 method for class 'corHerit'
\method{Initialize}{corHerit}(object, data, ...)
}
\arguments{
  \item{object}{
an object an object inheriting from class \code{corHerit} representing a correlation structure.
}
  \item{data}{a data frame in which to evaluate the variables defined in formula(object).
}
  \item{\dots}{this argument is included to make this method compatible with the generic.
}
}
\value{an initialized object with the same class as object representing a correlation structure.
}
\references{Pinheiro, J.C., and Bates, D.M. (2000) "Mixed-Effects Models in S and S-PLUS", Springer.
}

