\name{corHerit}
\alias{corHerit}
\title{Correlation Structure Constructor
}
\description{This function is a constructor for the \code{corHerit} class. Objects created using this constructor must later be initialized using the appropriate Initialize method.
}
\usage{
corHerit(value, paras, form, fixed)
}
\arguments{
  \item{value}{a vector with the initial values of the unknown parameters. 
}
  \item{paras}{a vector with additional known parameters.
}
  \item{form}{A one sided formula. When a
grouping factor is present in form, the correlation structure is assumed to apply only to observations within the same grouping level; observations with different grouping levels are assumed to be uncorrelated. Defaults to ~ 1, which corresponds to using the order of the observations in the data as a covariate, and no groups.

}
  \item{fixed}{an optional logical value indicating whether the coefficients should be allowed to vary in the optimization, or kept fixed at their initial value. Defaults to FALSE, in which case the coefficients are allowed to vary.
}
}
\value{an object of class \code{corHerit}, also inheriting from class \code{corStruct}, representing a correlation structure.
}
\references{
Pinheiro, J.C., and Bates, D.M. (2000) "Mixed-Effects Models in S and S-PLUS", Springer.
}



