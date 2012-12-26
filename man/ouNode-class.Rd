\name{ouNode-class}
\docType{class}
\alias{ouNode}
\alias{ouNode-class}
\title{Class \code{ouNode}}
\description{An object of class \code{"ouNode"} represents a motif node in a cluster tree}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ouNode", left, right, parent, distl, distr, sizel, sizer)}.
}
\section{Slots}{
  \describe{
    \item{\code{left}:}{character indicates the name of left leave}
    \item{\code{right}:}{character indicates the name of right leave}
    \item{\code{parent}:}{character indicates the name of parent node}
    \item{\code{distl}:}{numeric indicates the distance of left leave}
    \item{\code{distr}:}{numeric indicates the distance of right leave}
    \item{\code{sizel}:}{numeric indicates the size of left leave}
    \item{\code{sizer}:}{numeric indicates the size of right leave}
  }
}
\examples{
 new("ouNode", left="A", right="B", parent="Root", distl=1, distr=2, sizel=1, sizer=1)
}
\keyword{classes}