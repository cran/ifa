\name{ifa.bic}
\alias{ifa.bic}
\title{A function to compute the BIC}
\description{A function to compute the Bayesian Information Criterion (BIC), also known as Schwarz's Bayesian criterion (SBC), for the fitted IFA model, according to the formula \code{-2*log-likelihood + npar*log(nobs)}, where \code{npar} represents the number of parameters and \code{nobs} the number of observations in the fitted model.}
\usage{
ifa.bic(output)
}
\arguments{
  \item{output}{The fitted IFA model object, a list including the log-likelihood, the number of observations and the number of parameters.}
}
\value{
 It returns a numeric value with the corresponding BIC.
}


\references{Schwarz, G. (1978) Estimating the Dimension of a Model, Annals of Statistics, 6, 461-464.



Viroli, C. (2005). Choosing the number of factors in Independent Factor Analysis model, 
Metodoloski Zvezki, Advances in Methodology and Statistics, Vol. II, N. 2, 219-229.
Available at www2.stat.unibo.it/viroli.} 



\author{Cinzia Viroli}


\seealso{\link{ifa.aic}}
\examples{
data(memory)
fit<-ifa.em(memory$x,c(2,2),it=50,eps=0.001)
ifa.bic(fit)
}
\keyword{multivariate}
