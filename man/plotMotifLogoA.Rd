\name{plotMotifLogoA}
\alias{plotMotifLogoA}
\title{plot sequence logo without plot.new}
\description{
plot amino acid or DNA sequence logo in a given canvas
}
\usage{
plotMotifLogoA(pfm, font="Helvetica-Bold")
}
\arguments{
\item{pfm}{an object of pfm}
\item{font}{font of logo}
}
\value{
 none
}
\examples{
pcm<-matrix(runif(40,0,100),nrow=4,ncol=10)
pfm<-pcm2pfm(pcm)
rownames(pfm)<-c("A","C","G","T")
motif <- new("pfm", mat=pfm, name="bin_SOLEXA")
plotMotifLogoA(motif)
}