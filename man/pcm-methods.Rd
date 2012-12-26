\name{pcm-methods}
\docType{methods}
\alias{pcm}
\alias{addBlank,pcm,numeric,logical-method}
\alias{getIC,pcm,ANY-method}
\alias{pcm2pfm}
\alias{pcm2pfm,pcm,ANY-method}
\alias{pcm2pfm,matrix,ANY-method}
\alias{pcm2pfm,matrix,numeric-method}
\alias{pcm2pfm,data.frame,ANY-method}
\alias{pcm2pfm,data.frame,numeric-method}
\alias{matrixReverseComplement,pcm-method}
\alias{plot,pcm,ANY-method}

\title{"pcm" methods}
\description{
  methods for pcm objects. 
}

\usage{
	\S4method{addBlank}{pcm,numeric,logical}(x,n,b)
	\S4method{getIC}{pcm,ANY}(x,p="missing")
	\S4method{matrixReverseComplement}{pcm}(x)
	\S4method{plot}{pcm,ANY}(x,y="missing",\dots)
	\S4method{pcm2pfm}{pcm,ANY}(x,background="missing")
	\S4method{pcm2pfm}{matrix,ANY}(x,background="missing")
	\S4method{pcm2pfm}{matrix,numeric}(x,background)
	\S4method{pcm2pfm}{data.frame,ANY}(x,background="missing")
	\S4method{pcm2pfm}{data.frame,numeric}(x,background)
}

\arguments{
	\item{x}{An object of class \code{pfm}. For \code{getIC}, 
	if parameter p is followed, x should be an object of matrix.
	For \code{pcm2pfm}, x also could be an object of matrix.
	}
	\item{y}{Not use.}
	\item{p}{p is the background frequency.}
	\item{n}{how many spaces should be added.}
	\item{b}{logical value to indicate where the space should be added.}
	\item{background}{a \code{"numeric"} vector. The background frequency.}
	\item{\dots}{Further potential arguments passed to \code{plotMotifLogo}.}
}

\section{Methods}{
 \describe{
 	\item{addBlank}{\code{signature(x="pcm", n="numeric", b="logical")} add space 
  	into the position count matrix for alignment. b is a bool value, if TRUE, add 
  	space to the 3' end, else add space to the 5' end. n indicates how many 
  	spaces should be added.}
  	
    \item{coerce}{\code{signature(from = "pcm", to = "matrix")}: convert object pcm to matrix }
  	
  	\item{getIC}{\code{signature(x = "pcm",)} Calculate information
    content profile for position frequency matrix. }
    
    \item{matrixReverseComplement}{\code{signature(x = "pcm")} get the 
  	reverse complement of position frequency matrix.}
  	
	\item{plot}{\code{signature(x = "pcm")} Plots the
    sequence logo of the position count matrix. }
 } 
}
\examples{
pcm <- read.table(file.path(find.package("motifStack"), "extdata", "bin_SOLEXA.pcm"))
pcm <- pcm[,3:ncol(pcm)]
rownames(pcm) <- c("A","C","G","T")
motif <- new("pcm", mat=as.matrix(pcm), name="bin_SOLEXA")
getIC(motif)
matrixReverseComplement(motif)
as(motif,"matrix")
pcm2pfm(motif)
}
\keyword{classes}
