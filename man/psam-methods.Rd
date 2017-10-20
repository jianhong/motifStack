\name{psam-methods}
\docType{methods}
\alias{psam}
\alias{addBlank,psam,numeric,logical-method}
\alias{matrixReverseComplement,psam-method}
\alias{plot,psam,ANY-method}
\alias{$,psam-method}
\alias{$<-,psam-method}

\title{"psam" methods}
\description{
  methods for psam objects. 
}

\usage{
	\S4method{addBlank}{psam,numeric,logical}(x,n,b)
	\S4method{matrixReverseComplement}{psam}(x)
	\S4method{plot}{psam,ANY}(x,y="missing",\dots)
}

\arguments{
	\item{x}{An object of class \code{psam}.}
	\item{y}{Not use.}
	\item{n}{how many spaces should be added.}
	\item{b}{logical value to indicate where the space should be added.}
	\item{\dots}{Further potential arguments passed to \code{plotAffinityLogo}.}
}

\section{Methods}{
 \describe{
  \item{addBlank}{\code{signature(x="psam", n="numeric", b="logical")} add space 
  	into the position specific affinity matrix for alignment. b is a bool value, if TRUE, add 
  	space to the 3' end, else add space to the 5' end. n indicates how many 
  	spaces should be added.}

  \item{matrixReverseComplement}{\code{signature(x = "psam")} get the 
  	reverse complement of position specific affinity matrix.}
  	
  \item{plot}{\code{signature(x = "psam")} Plots the
    affinity logo of the position specific affinity matrix. }
  
  \item{$, $<-}{Get or set the slot of \code{\link{psam-class}}}
 } 
}
\examples{
motif <- importMatrix(file.path(find.package("motifStack"), "extdata", "PSAM.mxr"), 
                      format="psam")[[1]]
matrixReverseComplement(motif)
addBlank(motif, 1, FALSE)
addBlank(motif, 3, TRUE)
as(motif,"matrix")
}
\keyword{classes}
