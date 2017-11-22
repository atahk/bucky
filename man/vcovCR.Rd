\name{vcovCR}
\alias{vcovCR}
\title{
Clustered Robust Covariance Matrix Estimation
}
\description{
Robust estimation of the covariance matrix of the coefficient estimates
in regression models with clustering.
}
\usage{
vcovCR(x, cluster = NULL, type = c("CR", "CR0", "CR1"))
}
\arguments{
  \item{x}{
A fitted model object.
}
  \item{cluster}{
A variable or expression giving the cluster for each observation.
}
\item{type}{
A character string specifying the estimation type. For
details see below.
}
}
\details{
The default \code{type} of "CR" uses the same adjustment as Stata. The
values of "CR0" and "CR1" are analagous to "HC0" and "HC1",
respectively, in \code{vcovHC}.
}
\value{
A matrix containing the covariance matrix estimate with attribute
\code{type} giving the \code{type} option used in estimating it.
}

\seealso{
See Also \code{\link{robust.summary}} and \code{\link[sandwich]{vcovHC}}.
}
\examples{
clotting <- data.frame(
    cl = rep(1:2,each=9),
    u = c(5,10,15,20,30,40,60,80,100),
    lot = c(118,58,42,35,27,25,21,19,18,
            69,35,26,21,18,16,13,12,12))
clot.model <- glm(lot ~ log(u), data = clotting, family = Gamma)
vcovCR(clot.model, cluster=cl)

data(swiss)
model1 <- lm(Fertility ~ ., data = swiss)
## These should give the same answer
vcovCR(model1, cluster=1:nobs(model1), type="CR0")
vcovHC(model1, type="HC0")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ robust }% use one of  RShowDoc("KEYWORDS")
\keyword{ cluster }
