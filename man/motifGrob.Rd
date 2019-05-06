\name{motifGrob}
\alias{motifGrob}
\title{
Motif Grob
}
\description{
This function create a motif grob.
}
\usage{
motifGrob(pfm, x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                      width = unit(1, "npc"), height = unit(1, "npc"),
                      angle = 0, ic.scale=TRUE, default.units = "native", name=NULL, 
                      gp = gpar(fontfamily = "Helvetica-Bold",
                                fontface="bold"))
}
\arguments{
  \item{pfm}{an object of pfm}
  \item{x}{A numeric vector or unit object specifying x-values.}
  \item{y}{A numeric vector or unit object specifying y-values.}
  \item{width}{A numeric vector or unit object specifying width.}
  \item{height}{A numeric vector or unit object specifying height.}
  \item{angle}{A numeric value indicating the angle of rotation of the motif. Positive values indicate the amount of rotation, in degrees, anticlockwise from the positive x-axis.}
  \item{ic.scale}{logical If TRUE, the height of each column is proportional to its information content. Otherwise, all columns have the same height.}
  \item{default.units}{A string indicating the default units to use if x, y, width, or height are only given as numeric vectors.}
  \item{name}{A character value to uniquely identify the motifGrob once it has been pushed onto the grob tree.}
  \item{gp}{A gpar object, typically the output from a call to the function gpar. The list will be used as parameter of plotMotifLogoA.}
}
\value{
An gTree object.
}
\author{
Jianhong Ou
}
\examples{
pcm<-matrix(runif(40,0,100),nrow=4,ncol=10)
pfm<-pcm2pfm(pcm)
rownames(pfm)<-c("A","C","G","T")
motif <- new("pfm", mat=pfm, name="bin_SOLEXA")
motifGrob(motif)
}

