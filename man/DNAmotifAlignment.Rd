\name{DNAmotifAlignment}
\alias{DNAmotifAlignment}
\title{align DNA motifs}
\description{
align DNA motifs for plotting motifs stack
}
\usage{
DNAmotifAlignment(pfms, threshold=0.4, minimalConsensus=0, 
                  rcpostfix="(RC)", revcomp=rep(TRUE, length(pfms)))
}
\arguments{
\item{pfms}{a list of position frequency matrices, pfms must be a list of class pfm}
\item{threshold}{information content cutoff threshold for useful postions}
\item{minimalConsensus}{minimal length of consensus for alignment}
\item{rcpostfix}{the postfix for reverse complements}
\item{revcomp}{a logical vector to indicates whether the reverse complemet should be involved into alignment}
}
\value{
 a list of aligned motifs
}
\examples{
pcms<-readPCM(file.path(find.package("motifStack"), "extdata"),"pcm$")
motifs<-lapply(pcms,pcm2pfm)
motifs<-DNAmotifAlignment(motifs)
}