\name{plotMotifLogo}
\alias{plotMotifLogo}
\title{plot sequence logo}
\description{
plot amino acid or DNA sequence logo
}
\usage{
plotMotifLogo(pfm, motifName, p=rep(0.25, 4), font="Helvetica-Bold", 
colset=c("#00811B","#2000C7","#FFB32C","#D00001"), 
xaxis=TRUE,yaxis=TRUE,xlab="position",ylab="bits",
xlcex=1.2, ylcex=1.2, ncex=1.2, ic.scale=TRUE)
}
\arguments{
\item{pfm}{a position frequency matrices}
\item{motifName}{motif name}
\item{p}{background possibility}
\item{font}{font of logo}
\item{colset}{color setting for each logo letter}
\item{xaxis}{draw x-axis or not}
\item{yaxis}{draw y-axis or not}
\item{xlab}{x-label, do nothing if set xlab as NA}
\item{ylab}{y-label, do nothing if set ylab as NA}
\item{xlcex}{cex value for x-label}
\item{ylcex}{cex value for y-label}
\item{ncex}{cex value for motif name}
\item{ic.scale}{logical If TRUE, the height of each column is proportional to its information content. Otherwise, all columns have the same height.}
}
\value{
 none
}
\examples{
pcm<-matrix(runif(40,0,100),nrow=4,ncol=10)
pfm<-pcm2pfm(pcm)
rownames(pfm)<-c("A","C","G","T")
plotMotifLogo(pfm)
}