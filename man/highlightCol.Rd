\name{highlightCol}
\alias{highlightCol}
\title{
add alpha transparency value to a color
}
\description{
An alpha transparency value can be specified to a color, in order to get better
color for background.
}
\usage{
highlightCol(col, alpha = 0.5)
}
\arguments{
  \item{col}{
    vector of any of the three kinds of R color specifications, 
    i.e., either a color name (as listed by \link[grDevices]{colors}()), 
    a hexadecimal string of the form "#rrggbb" or "#rrggbbaa" (see \link[grDevices]{rgb}), 
    or a positive integer i meaning \link[grDevices]{palette}()[i].}
  \item{alpha}{a value in [0, 1]}
}
\value{
a vector of colors in hexadecimal string of the form "#rrggbbaa".
}
\author{
Jianhong Ou
}
\examples{
    highlightCol(1:5, 0.3)
    highlightCol(c("red", "green", "blue"), 0.3)
}
\keyword{ misc }
