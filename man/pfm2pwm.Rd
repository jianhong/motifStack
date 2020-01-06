\name{pfm2pwm}
\alias{pfm2pwm}
\title{
convert pfm object to PWM
}
\description{
convert pfm object to PWM
}
\usage{
pfm2pwm(x, N=10000)
}
\arguments{
  \item{x}{
 an object of \code{\link{pfm}} or \code{\link{pcm}} or matrix
}
  \item{N}{
  Total number of event counts used for pfm generation.
  }
}
\value{
A numeric matrix representing the Position Weight Matrix for PWM.
}
\seealso{
\code{\link[Biostrings:matchPWM]{PWM}}
}
\author{
 Jianhong Ou
}
\examples{
library("MotifDb")
matrix.fly <- query(MotifDb, "Dmelanogaster")
pfm2pwm(matrix.fly[[1]])
}
\keyword{ misc }