\name{psam-class}
\docType{class}
\alias{psam-class}

\title{Class \code{"psam"}}
\description{
  An object of class \code{"psam"} represents the position specific affinity
  matrix (PSAM) of a DNA/RNA/amino-acid sequence motif. The entry stores a matrix,
  which in row i, column j gives the
  affinity of observing nucleotide/or amino acid i in
  position j of the motif.  
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("psam", mat, name, alphabet, color)}.
}
\section{Slots}{
  \describe{
    \item{\code{mat}}{Object of class \code{"matrix"} The position specific affinity matrix}
    \item{\code{name}}{code{"character"} The motif name}
    \item{\code{alphabet}}{\code{"character"} The sequence
      alphabet. "DNA", "RNA", "AA" or "others".}
    \item{\code{color}}{a \code{"character"} vector. The color setting for each symbol}
    \item{\code{tags}}{a list to save metadata of motif.}
    \item{\code{markers}}{a list of marker object.}
  }
}
\section{Methods}{
 \describe{
	\item{addBlank}{\code{signature(x="psam", n="numeric", b="logical")} add space 
  	into the position specific affinity matrix for alignment. b is a bool value, if TRUE, add 
  	space to the 3' end, else add space to the 5' end. n indicates how many 
  	spaces should be added.}
  	
  	\item{coerce}{\code{signature(from = "psam", to = "matrix")}: convert object psam to matrix }
    
    \item{matrixReverseComplement}{\code{signature(x = "psam")} get the 
  	reverse complement of position specific affinity matrix.}
  	
	\item{plot}{\code{signature(x = "psam")} Plots the
    affinity logo of the position specific affinity matrix. }
 } 
}
\examples{
motif <- importMatrix(file.path(find.package("motifStack"), "extdata", "PSAM.mxr"), 
         format="psam")[[1]]
plot(motif)
}
\keyword{classes}
