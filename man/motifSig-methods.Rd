\name{motifSig-methods}
\docType{methods}
\alias{motifSig}
\alias{signatures}
\alias{signatures,motifSig-method}
\alias{frequence}
\alias{frequence,motifSig-method}
\alias{nodelist}
\alias{nodelist,motifSig-method}

\title{"motifSig" methods}
\description{
  methods for motifSig objects. 
}

\usage{
	\S4method{signatures}{motifSig}(object)
	\S4method{frequence}{motifSig}(object)
	\S4method{nodelist}{motifSig}(object)
}

\arguments{
	\item{object}{An object of class \code{motifSig}.}
}
\section{Methods}{
 \describe{
   	\item{signatures}{\code{signature(object = "motifSig")} return the signatures of motifSig }
    
    \item{frequence}{\code{signature(object = "motifSig")} return the frequency of motifSig}
  	
	\item{nodelist}{\code{signature(object = "motifSig")} return the nodelist of motifSig}
 } 
}
\keyword{classes}
