\name{mergeMotifs}
\alias{mergeMotifs}
\title{
    merge multiple motifs
}
\description{
    merge multiple motifs by calculate mean of each position
}
\usage{
mergeMotifs(..., bgNoise=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{\link{pcm} or \link{pfm} objects}
  \item{bgNoise}{if it is not NA, test will using a background by Dirichlet(1)-distributed random frequencies with weight bg.noise. 
  The value of bgNoise should be a number in the range of 0 to 1, eg. 0.05}
}
\value{
    a \link{pfm} object
}
\author{
Jianhong Ou
}
\examples{
    pcms<-readPCM(file.path(find.package("motifStack"), "extdata"),"pcm$")
    mergeMotifs(pcms)
}
\keyword{ misc }
