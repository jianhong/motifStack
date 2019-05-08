\name{marker-class}
\docType{class}
\alias{marker}
\alias{marker-class}
\alias{$,marker-method}
\alias{$<-,marker-method}
\title{Class \code{marker}}
\description{An object of class \code{"marker"} represents a marker in a motif}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("marker", type, start, stop, label, gp)}.
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{character(1), 1 of 'rect', 'text', 'line'}
    \item{\code{start}:}{numberic indicates the start position of the marker}
    \item{\code{stop}:}{numberic indicates the end position of the marker}
    \item{\code{label}:}{character(1) indicates the text to be marked}
    \item{\code{gp}:}{an object of gpar}
  }
}
\section{Methods}{
    \describe{
        \item{$, $<-}{Get or set the slot of \code{\link{marker}}}
    }
}
\examples{
 new("marker", type="rect", start=c(2, 4), gp=gpar(lty=3))
}
\keyword{classes}