\name{pfm-class}
\docType{class}
\alias{pfm-class}

\title{Class \code{"pfm"}}
\description{
  An object of class \code{"pfm"} represents the position frequency
  matrix of a DNA/RNA/amino-acid sequence motif. The entry stores a matrix,
  which in row i, column j gives the
  frequency of observing nucleotide/or amino acid i in
  position j of the motif.  
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("pfm", mat, name, alphabet, color, background)}.
}
\section{Slots}{
  \describe{
    \item{\code{mat}}{Object of class \code{"matrix"} The position frequency matrix}
    \item{\code{name}}{code{"character"} The motif name}
    \item{\code{alphabet}}{\code{"character"} The sequence
      alphabet. "DNA", "RNA", "AA" or "others".}
    \item{\code{color}}{a \code{"character"} vector. The color setting for each symbol}
    \item{\code{background}}{a \code{"numeric"} vector. The background frequency.}
  }
}
\section{Methods}{
 \describe{
	\item{addBlank}{\code{signature(x="pfm", n="numeric", b="logical")} add space 
  	into the position frequency matrix for alignment. b is a bool value, if TRUE, add 
  	space to the 3' end, else add space to the 5' end. n indicates how many 
  	spaces should be added.}
  	
  	\item{coerce}{\code{signature(from = "pfm", to = "matrix")}: convert object pfm to matrix }
  	
  	\item{getIC}{\code{signature(x = "pfm",)} Calculate information
    content profile for position frequency matrix. }

	\item{getIC}{\code{signature(x = "matrix", p = "numeric")} Calculate
    information content profile for matrix. p is the background frequency}
    
    \item{matrixReverseComplement}{\code{signature(x = "pfm")} get the 
  	reverse complement of position frequency matrix.}
  	
  	\item{trimMotif}{\code{signature(x = "pfm", t= "numeric")} trim motif by information content. }
  	
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
plot(motif)
}
\keyword{classes}
