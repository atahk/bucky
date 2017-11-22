\name{summary.mi.estimates}
\alias{summary.mi.estimates}
\title{
Summary for multiple imputation
}
\description{
Output summary information for estimates computed on multiply imputed data sets.
}
\usage{
\method{summary}{mi.estimates}(object, \dots)
}
\arguments{
  \item{object}{
An object of class \code{mi.estimates}, usually computed with
\code{\link{mi.eval}}.
}
  \item{\dots}{
Any additional arguments to be passed to \code{coeftest}.
}
}
\details{
  This function provides summary output for models estimated
  using multiple imputation using \code{\link{mi.eval}}.
}
\value{
An object of class \code{summary.mi.estimates} containing a
\code{coefficients} object of class \code{coeftest} and other
summary information.
}

\seealso{
See Also \code{\link{mi.eval}} and \code{\link[Amelia]{amelia}}.
}
\examples{
if (require("Amelia")) {
    data(africa)
    a.out <- amelia(x = africa, cs = "country", ts = "year", logs = "gdp_pc")

    model <- mi.eval(lm(civlib ~ log(gdp_pc), data=a.out, subset=year==1973))

    summary(model)
}
}
\keyword{ model }% use one of  RShowDoc("KEYWORDS")
\keyword{ htest }
