\name{pfm-methods}
\docType{methods}
\alias{pfm}
\alias{addBlank}
\alias{addBlank,pfm,numeric,logical-method}
\alias{getIC}
\alias{getIC,pfm,ANY-method}
\alias{getIC,matrix,numeric-method}
\alias{matrixReverseComplement}
\alias{matrixReverseComplement,pfm-method}
\alias{plot,pfm-method}

\title{"pfm" methods}
\description{
  methods for pfm objects. 
}

\usage{
	\S4method{addBlank}{pfm,numeric,logical}(x,n,b)
	\S4method{getIC}{pfm,ANY}(x,p="missing")
	\S4method{getIC}{matrix,numeric}(x,p)
	\S4method{matrixReverseComplement}{pfm}(x)
	\S4method{plot}{pfm}(x,y="missing",\dots)
}

\arguments{
	\item{x}{An object of class \code{pfm}. For \code{getIC}, 
	if parameter p is followed, x should be an object of matrix.}
	\item{y}{Not use.}
	\item{p}{p is the background frequency.}
	\item{n}{how many spaces should be added.}
	\item{b}{logical value to indicate where the space should be added.}
	\item{\dots}{Further potential arguments passed to \code{plotMotifLogo}.}
}

\section{Methods}{
 \describe{
  \item{addBlank}{\code{signature(x="pfm", n="numeric", b="logical")} add space 
  	into the position frequency matrix for alignment. b is a bool value, if TRUE, add 
  	space to the 3' end, else add space to the 5' end. n indicates how many 
  	spaces should be added.}
  
  \item{getIC}{\code{signature(x = "pfm",)} Calculate information
    content profile for position frequency matrix. }

  \item{getIC}{\code{signature(x = "matrix", p = "numeric")} Calculate
    information content profile for matrix. p is the background frequency}

  \item{matrixReverseComplement}{\code{signature(x = "pfm")} get the 
  	reverse complement of position frequency matrix.}
  	
  \item{plot}{\code{signature(x = "pfm")} Plots the
    sequence logo of the position frequency matrix. }
 } 
}
\examples{
pcm <- read.table(file.path(find.package("motifStack"), "extdata", "bin_SOLEXA.pcm"))
pcm <- pcm[,3:ncol(pcm)]
rownames(pcm) <- c("A","C","G","T")
motif <- pcm2pfm(pcm)
motif <- new("pfm", mat=motif, name="bin_SOLEXA")
getIC(motif)
matrixReverseComplement(motif)
addBlank(motif, 1, FALSE)
addBlank(motif, 3, TRUE)
as(motif,"matrix")
}
\keyword{classes}
