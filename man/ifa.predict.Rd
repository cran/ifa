\name{ifa.predict}
\alias{ifa.predict}
\title{A function to predict the latent variables}
\description{
 A function to compute the reconstructed factors.
}
\usage{
ifa.predict(y, output, method = "lms")
}
\arguments{
  \item{y}{The data matrix for which factor reconstructions are desired}
  \item{output}{The result of a call to ifa.em}
  \item{method}{The method for reconstructing the factors}
}
\details{
 The function reconstructs the factor scores by different methods: \code{"bartlett"}, \code{"thompson"} 
 or \code{"lms"} (least mean squares). The default is \code{"lms"}.
}
\value{
  A matrix with factors in columns.
}

\author{Cinzia Viroli}
\keyword{multivariate}

