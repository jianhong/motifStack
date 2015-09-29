\name{plotMotifLogoA}
\alias{plotMotifLogoA}
\title{plot sequence logo without plot.new}
\description{
plot amino acid or DNA sequence logo in a given canvas
}
\usage{
plotMotifLogoA(pfm, font="Helvetica-Bold", ic.scale=TRUE, fontsize=12)
}
\arguments{
\item{pfm}{an object of pfm}
\item{font}{font of logo}
\item{ic.scale}{logical If TRUE, the height of each column is proportional to its information content. Otherwise, all columns have the same height.}
\item{fontsize}{font size of the template for grImport, default 12. Higher value make better quality figure, but also increase the file size.}
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