\name{venus}
\alias{venus}
\alias{plot.venus}
\alias{print.venus}
\alias{summary.venus}
\alias{coef.venus}
\alias{effects.venus}
\alias{fitted.venus}
\alias{print.summary.venus}
\alias{residuals.venus}
\alias{vcov.venus}
\title{
Do a two-stage multivariate adaptive regression splines (MARS) fit.
}
\description{
When some parameters are nuisance parameters, it may be desirable to regress them out of the model first.  This
function does this in a combination of multiple regression
and multivariate adaptive regression splines.
}
\usage{
venus(formula, nuisance = NULL, data = NULL,
      method = c("earth", "lm"),
      nuisanceArgs = NULL, mainArgs = NULL, 
      svd.thresh=.99, ps_method = c("po", "ds"), ...)
\method{plot}{venus}(x, types = c("nuisance", "main"), ...)
\method{print}{venus}(x, types = c("nuisance", "main"), ...)
\method{summary}{venus}(object, ...)
\method{coef}{venus}(object, type = "main", ...)
\method{effects}{venus}(object, type = "main", ...)
\method{fitted}{venus}(object, type = "main", ...)
\method{print}{summary.venus}(x, types = c("nuisance", "main"), ...)
}
\arguments{
  \item{formula}{
The main model containing the parameters of interest.
}
  \item{nuisance}{
The nuisance model; predictors here need to be regressed out.
}
  \item{data}{
A dataframe to fit.
}
  \item{method}{
A vector with two entries from \code{c("earth", "lm")}
giving the methods to be used in the two stages.
}
  \item{nuisanceArgs}{
A list of additional arguments to pass to the nuisance fit.
}
  \item{mainArgs}{
A list of additional arguments to pass to the main fit.
}
  \item{x, object}{ Object for S3 methods.}
  \item{type, types}{ Which part of the \code{"venus"}
  object should be used?  \code{type} allows one of
  \code{c("main", "nuisance", "predictor")}, \code{types}
  allows more than one from the same list.}
  \item{svd.thresh}{Threshold on the amount of variance captured in the singular value decomposition of the control variable matrix. }
  \item{ps_method}{Method for controlling for nuisance variables can be one of \code{"po"} for partialling out or \code{"ds"} for double selection.}
  \item{...}{Additional parameters to pass to component methods.}
}
\details{
We use the \code{\link{earth}} function from the package
of the same name to do the MARS fits.  This may be done
at either the first or second stage (or both), depending on
the choice of values in \code{method}; the default is
MARS in the first stage, linear fit in the second stage.

The algorithm is as follows:

\enumerate{
\item The nuisance model is fit using \code{method[1]} (\code{\link{earth}} by default).  This results in the
choice of predictors after adaptation to the data.
\item The chosen model is used to find residuals in the
response and in the predictors in the main model.
\item The residuals in the response are fit against the
residuals in the main predictors using \code{method[2]}
(\code{\link{lm}} by default).
\item Both fitted models are returned.
}

}
\value{
\code{venus} returns a list with S3 class \code{"venus"} containing three elements:
\item{nuisanceFit}{The first stage fit}
\item{mainFit}{The second stage fit}
\item{predictorFit}{The linear fit to partial the nuisance
variables out of the predictors in the main fit}

The common methods for linear models are defined for these 
objects; each has a \code{type} or \code{types} argument to
select which of the three components is/are used in the call.
For example, \code{coef(venusObject, type = "nuisance")} would extract the
coefficients from the nuisance fit.
}
\author{
Duncan Murdoch and David Armstrong
}

\examples{
  if (!require(MASS))
    stop("This example requires MASS")

  sigma2.e.x <- 4.7*3
  ## sigma2.e is the residual variance for y. On average, the systematic
  ## variance in y owing to x and Z is around 7.3.  Initially, I make the
  ## residual variance half of that figure.
  sigma2.e <- 7.3/2
  ## Z are the control variables or nuisance variables
  Z <- mvrnorm(1000, c(0,0,0),
             matrix(c(1,.3,.3,.3,1,.3,.3,.3,1), ncol=3),
             empirical=TRUE)
  ## X is the variable of interest that has a non-linear relationship to
  ## one of the elements in Z (though this isn't strictly speaking
  ## necessary, it does highlight one of the nice features of the
  ## MARS approach, which is that you don't need to know the functional
  ## form of the relationship of controls to variables of interest
  ## (e.g., y and X)
  X <- Z[,1] + Z[,1]^2 + Z[,2] + rnorm(1000, 0, sqrt(sigma2.e.x))
  ## collect results in a data frame
  df <- data.frame(z1 = Z[,1], z2=Z[,2], z3 = Z[,3], x=scale(X))
  ## make y a linear function of x and the z variables,
  ## so the theoretical coefficient on x is 1
  df$y <- with(df, z1 + z2 + z3 + x) + rnorm(1000, 0, sqrt(sigma2.e))

  venus(y ~ x, nuisance = y ~ .-x, data = df)

  # Compare to the pure linear fit

  venus(y ~ x, nuisance = y ~ .-x, method = c("lm", "lm"), data=df)
  lm(y ~ ., data = df)
}
