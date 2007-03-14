\name{ifa.init.random}
\alias{ifa.init.random}
\title{A function to randomly inizialize the IFA parameters}
\description{
  A function to initialize the IFA parameters to random values.
}
\usage{
ifa.init.random(y,L)
}
\arguments{
  \item{y}{The data to be fitted}
  \item{L}{The number of factors}
}
\value{A list containing 2 initialized parameters:
  \item{H}{The factor loading matrix}
  \item{psi}{The noise diagonal variance matrix}
}

\seealso{\code{\link{ifa.init.pca}}}
\examples{
data(memory)
init.values<-ifa.init.random(memory$x,2)
fit<-ifa.em(memory$x,c(2,2),it=50,eps=0.0001,init.values)
}

\keyword{multivariate}
