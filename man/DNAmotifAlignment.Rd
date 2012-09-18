\name{DNAmotifAlignment}
\alias{DNAmotifAlignment}
\title{align DNA motifs}
\description{
align DNA motifs for plotting motifs stack
}
\usage{
DNAmotifAlignment(pfms, threshold=0.4)
}
\arguments{
\item{pfms}{a list of position frequency matrices, pfms must be a list of class pfm}
\item{threshold}{information content cutoff threshold for useful postions}
}
\value{
 a list of aligned motifs
}
\examples{
pcms<-readPCM(file.path(find.package("motifStack"), "extdata"),"pcm$")
pcms<-lapply(pcms,function(.ele){.ele<-.ele[,3:ncol(.ele)];rownames(.ele)<-c("A","C","G","T");.ele})
motifs<-lapply(pcms,pcm2pfm)
motifs<-lapply(names(motifs), function(.ele, motifs){new("pfm",mat=motifs[[.ele]], name=.ele)},motifs)
motifs<-DNAmotifAlignment(motifs)
}