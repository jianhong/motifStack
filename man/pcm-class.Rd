\name{pcm-class}
\docType{class}
\alias{pcm-class}

\title{Class \code{"pcm"}}
\description{
  An object of class \code{"pcm"} represents the position count
  matrix of a DNA/RNA/amino-acid sequence motif. The entry stores a matrix,
  which in row i, column j gives the
  counts of observing nucleotide/or amino acid i in
  position j of the motif.  
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("pcm", mat, name, alphabet, color, background)}.
}
\section{Slots}{
  \describe{
    \item{\code{mat}}{Object of class \code{"matrix"} The position count matrix}
    \item{\code{name}}{code{"character"} The motif name}
    \item{\code{alphabet}}{\code{"character"} The sequence
      alphabet. "DNA", "RNA", "AA" or "others".}
    \item{\code{color}}{a \code{"character"} vector. The color setting for each symbol}
    \item{\code{background}}{a \code{"numeric"} vector. The background frequency.}
  }
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
  	
  	\item{trimMotif}{\code{signature(x = "pcm", t= "numeric")} trim motif by information content. }
  	
	\item{plot}{\code{signature(x = "pcm")} Plots the
    sequence logo of the position count matrix. }
 } 
}
\examples{
pcm <- read.table(file.path(find.package("motifStack"), "extdata", "bin_SOLEXA.pcm"))
pcm <- pcm[,3:ncol(pcm)]
rownames(pcm) <- c("A","C","G","T")
motif <- new("pcm", mat=as.matrix(pcm), name="bin_SOLEXA")
plot(motif)
}
\keyword{classes}
