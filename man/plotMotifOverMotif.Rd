\name{plotMotifOverMotif}
\alias{plotMotifOverMotif}
\title{plot motif over another motif}
\description{
plot motif over another motif to emphesize the difference.
}
\usage{
plotMotifOverMotif(motif, backgroundMotif, bgNoise=NA, 
                font="Helvetica-Bold", textgp=gpar())
}
\arguments{
\item{pfm}{an object of \code{\link{pcm}} or \code{\link{pfm}}}
\item{backgroundMotif}{an object of \code{\link{pcm}} or \code{\link{pfm}}}
\item{bgNoise}{if it is not NA, test will using a background by Dirichlet(1)-distributed random frequencies with weight bg.noise. 
  The value of bgNoise should be a number in the range of 0 to 1, eg. 0.05}
\item{font}{font for logo symbol}
\item{textgp}{text parameter}
}
\value{
 none
}
\examples{
pcms <- readPCM(file.path(find.package("motifStack"), "extdata"),"pcm$")
len <- sapply(pcms, function(.ele) ncol(.ele$mat))
pcms <- pcms[len==7]
plotMotifOverMotif(pcms[[1]], pcms[[2]], bgNoise=0.05)
}