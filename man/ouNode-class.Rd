% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signature.R
\docType{class}
\name{ouNode-class}
\alias{ouNode-class}
\alias{ouNode}
\alias{$,ouNode-method}
\alias{$<-,ouNode-method}
\title{Class \code{ouNode}}
\usage{
\S4method{$}{ouNode}(x, name)
}
\arguments{
\item{x}{A ouNode object}

\item{name}{slot name of ouNode object}
}
\description{
An object of class \code{"ouNode"} represents a motif node in a cluster tree
}
\section{Objects from the Class}{
 Objects can be created by calls of the form
\code{new("ouNode", left, right, parent, distl, distr, sizel, sizer)}.
}

\examples{

 new("ouNode", left="A", right="B", parent="Root", distl=1, distr=2, sizel=1, sizer=1)

}
\keyword{classes}
