\name{pcm2pfm}
\alias{pcm2pfm}
\title{convert pcm to pfm}
\description{
convert position count matrix to position frequency matrix
}
\usage{
pcm2pfm(pcm)
}
\arguments{
\item{pcm}{position count matrix}
}
\value{
 position weight matrix
}
\examples{
pcms<-readPCM(file.path(find.package("motifStack"), "extdata"),"pcm$")
pcms<-lapply(pcms,function(.ele){.ele<-.ele[,3:ncol(.ele)];rownames(.ele)<-c("A","C","G","T");.ele})
motifs<-lapply(pcms,pcm2pfm)

}