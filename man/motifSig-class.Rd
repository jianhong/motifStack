\name{motifSig-class}
\docType{class}
\alias{motifSig-class}

\title{Class \code{"motifSig"}}
\description{
  An object of class \code{"motifSig"} represents the output of function \link{motifSignature}
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("motifSig", signature, freq, nodelist, gpcol)}.
}
\section{Slots}{
  \describe{
    \item{\code{signatures}}{list object of class \code{"pfm"}}
    \item{\code{freq}}{code{"numeric"} signature frequency}
    \item{\code{nodelist}}{list object of class \code{"\link{ouNode}"}}
    \item{\code{gpcol}}{code{"character"} signature group color sets}
  }
}
\section{Methods}{
 \describe{
   	\item{signatures}{\code{signature(object = "motifSig")} return the signatures of motifSig }
    
    \item{frequence}{\code{signature(object = "motifSig")} return the frequency of motifSig}
  	
	\item{nodelist}{\code{signature(object = "motifSig")} return the nodelist of motifSig}
    
    \item{sigColor}{\code{signature(object = "motifSig")} return the group color sets of motifSig}
 } 
}
\keyword{classes}
