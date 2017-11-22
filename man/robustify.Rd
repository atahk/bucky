\name{robustify}
\alias{robustify}
\title{
Robustify a model
}
\description{
Create a robustified object that includes robust or clustered robust
standard errors.
}
\usage{
robustify(x, cluster, type, omega, ...)
}
\arguments{
  \item{x}{
A model of class \code{lm}, \code{glm}, or any other class which
contains a \code{call} object and methods for \code{estfun} and \code{nobs}.
}
  \item{cluster}{
The variable on which to cluster (if any). If this is not specified,
unclustered robust standard errors using \code{\link[sandwich]{vcovHC}} are
used. If this is specified, clustered robust standard errors using
\code{\link{vcovCR}} are used.
}
  \item{type}{
A character string specifying the estimation type. The default is to use the
defaults for \code{\link[sandwich]{vcovHC}} or \code{\link{vcovCR}}, depending on whether
\code{cluster} is specified. For details, see \code{\link[sandwich]{vcovHC}} or \code{\link{vcovCR}}.
}
  \item{omega}{
A vector or a function depending on the arguments `residuals'
(the working residuals of the model), `diaghat' (the diagonal
of the corresponding hat matrix) and `df' (the residual
degrees of freedom). For details, see \code{\link[sandwich]{vcovHC}} or \code{\link{vcovCR}}.
}
  \item{\dots}{
Any additional arguments to be passed to \code{\link[lmtest]{coeftest}}.
}
}
\details{
This function creates a \code{robustified} object containing the
model, coefficients, and variance-covariance matrix based on
\code{\link[sandwich]{vcovHC}} or \code{\link{vcovCR}}, respectively. If the \code{cluster}
option is specified, the variance-covariance matrix is computed using
\code{\link{vcovCR}} with clustering on \code{cluster}. If not, the
variance-covariance matrix is computed using
\code{\link[sandwich]{vcovHC}}.
}
\value{
An object of class \code{robustified} with the \code{method}
attribute specifying the type of standard errors used.
}

\seealso{
See Also \code{\link{summary.robustified}}, \code{\link[sandwich]{vcovHC}}, \code{\link{vcovCR}} and \code{\link[lmtest]{coeftest}}.
}
\examples{
## With clustering
clotting <- data.frame(
    cl = 1:9,
    u = c(5,10,15,20,30,40,60,80,100),
    lot = c(118,58,42,35,27,25,21,19,18,
            69,35,26,21,18,16,13,12,12))
clot.model <- glm(lot ~ log(u), data = clotting, family = Gamma)
robust.clot.model <- robustify(clot.model, cluster=cl)
coeftest(robust.clot.model)

## Without clustering
data(swiss)
model1 <- robustify(lm(Fertility ~ ., data = swiss))
coeftest(model1)
}
\keyword{ robust }% use one of  RShowDoc("KEYWORDS")
\keyword{ cluster }
\keyword{ htest }
