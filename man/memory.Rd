\name{memory}
\alias{memory}
\docType{data}
\title{Memory Tests on Alzheimer and Depressed Patients}
\description{
 Data from study by Robert Hart, concerning ability of memory tests
to distinguish between patients with mild Alzheimer's disease,
depression, and controls.
}
\usage{data(memory)}

\format{memory is a list of 2 components:
 \code{$y} is the group variable (1=Alzheimer; 2=Depressed; 3=Normal) and \code{$x} is a data.frame containing 45 obs. 
 of 6 variables: 
 
 HIR : HI Imagery Recall 
 
 LIR : LO Imagery Recall
 
 HIUM: HI Imagery Unreminded  Memory
  
 LIUM: LO Imagery Unreminded Memory 
 
 SM  : Storage 
 
 RM  : Recognition Memory 
}

\references{
  Hart, R.P., Kwentus, J.A., Taylor, J.R., \& Hamer, R.M. (1988), Productive naming and
memory in depression and Alzheimer's type dementia, Archives
of Clinical Neuropsychology, 3, 313-322.
}
\examples{
data(memory)
}
\keyword{datasets}
