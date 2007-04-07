\name{ifa.aic}
\alias{ifa.aic}
\title{A function to compute the AIC}
\description{
 A function to compute the Akaike Information Criterion (AIC) for the fitted IFA model, according to the formula \code{-2*log-likelihood + 2*npar}, where \code{npar} represents the number of parameters.}
\usage{
ifa.aic(output)
}
\arguments{
  \item{output}{The fitted IFA model object, a list including the log-likelihood and the number of parameters}
}
\value{
  It returns a numeric value with the corresponding AIC.
}


\references{Sakamoto, Y., Ishiguro, M., and Kitagawa G. (1986). Akaike Information Criterion Statistics. D. Reidel Publishing Company.


Viroli, C. (2005). Choosing the number of factors in Independent Factor Analysis model, 
Metodoloski Zvezki, Advances in Methodology and Statistics, Vol. II, N. 2, 219-229.
Available at \url{www2.stat.unibo.it/viroli}.} 


\author{Cinzia Viroli}


\seealso{\link{ifa.bic}}
\examples{
data(memory)
fit<-ifa.em(memory$x,c(2,2),it=50,eps=0.001)
ifa.aic(fit)
}
\keyword{multivariate}


