\name{ifa.em}
\alias{ifa.em}
\title{Fitting an Independent Factor Analysis model by the EM algorithm}
\description{
  \code{ifa.em} is used to perform Independent Factor Analysis on a matrix of data by the Expectation Maximization algorithm.  
}
\usage{
ifa.em(y, ni, it = 15, eps = 0.001, init = NULL, scaling = TRUE)
}
\arguments{
  \item{y}{A data matrix with \code{n} rows representing observations and \code{p} columns representing variables.}
  \item{ni}{A vector containing the number of mixture components for modeling each factor. The number of factors is equal to the length of this vector.}
  \item{it}{The maximum number of iterations of the EM algorithm. By default it is set to 15.}
  \item{eps}{The lower bound for relative variation of the likelihood. It is used as alternative stopping rule for the EM algorithm: if the relative increment of the likelihood is lower than \code{eps} the fitting is stopped. The default is 0.001.}
  \item{init}{A list containing initial values for the factor loading matrix (specified by \code{H}) and the variance matrix of the noise term (specified by \code{psi}). If \code{NULL} (default) the algorithm starts from the ordinal factor model solution.}
  \item{scaling}{If TRUE (default) the data are scaled before fitting the IFA model}
}
\details{
Independent Factor Analysis is a latent variable model with
independent and non gaussian factors. The \code{p} observed
variables \code{x} are modelled in terms of a smaller set of \code{k}
unobserved independent latent variables, \code{y}, and an additive
specific term u: x=Hy+u,
where \code{u} is assumed
to be normally distributed with diagonal variance matrix \code{Psi} and the factor loading
matrix \code{H} is also
termed as \code{mixing matrix}. The density of each 
factor is modelled by a mixture of gaussians. The model is fitted by the EM algorithm. 
This version can be computationally slow in the complex cases. 
A faster R package for window which is based on fortran code can be downloaded at the home: 
www2.stat.unibo.it/viroli in the section Software and Data.}
\value{A list containing the components:
  \item{H}{The estimated factor loading matrix}
  \item{lik}{The log-likelihood computed at each iteration of the EM algorithm}
  \item{w}{A matrix with the estimated weights of the mixtures. Each row contains the weights of each factor}
  \item{mu}{A matrix with the estimated component means of the mixtures. Each row contains the vector means for each factor}
  \item{vu}{A matrix with the estimated component variances of the mixtures. Each row contains the vector variances for each factor}
  \item{psi}{The noise diagonal variance matrix}
  \item{ni}{The input vector specifying the number of components for each mixture}
  \item{L}{The number of factors}
  \item{numvar}{The number of observed variables}
  \item{numobs}{The number of observations}
}
\references{
Attias H. (1999), Independent Factor Analysis, Neural Computation, 11, 803--851.

Montanari A. and Viroli C. (2010), The Independent Factor Analysis 
approach to latent variable modeling, Statistics, 44, 397--416. 
}
\author{Cinzia Viroli}

\seealso{\code{\link{ifa.init.random}}, \code{\link{ifa.init.pca}}}
\examples{
data(memory)
init.values<-ifa.init.random(memory$x,2)
fit<-ifa.em(memory$x,c(2,2),it=50,eps=0.0001,init.values)

fit<-ifa.em(memory$x,c(2,2),it=50,eps=0.0001)
plot(ifa.predict(scale(memory$x),fit))
}

\keyword{multivariate}
