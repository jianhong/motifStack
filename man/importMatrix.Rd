\name{importMatrix}
\alias{importMatrix}
\title{
import motifs from local files
}
\description{
Import the motifs into \link{pcm-class} or \link{pfm-class} from files exported from 
Transfac, CisBP, and JASPAR. 
}
\usage{
importMatrix(filenames, 
             format=c("auto", "pfm", "cm", "pcm", "meme", 
                      "transfac", "jaspar", "scpd", "beeml", "cisbp",
                      "psam"), 
             to=c("auto", "pcm", "pfm", "psam"))
}
\arguments{
  \item{filenames}{
filename to be imported.
}
  \item{format}{
file format
}
  \item{to}{
import to \link{pcm-class} or \link{pfm-class}
}
}
\value{
a list of object \link{pcm-class} or \link{pfm-class}
}
\author{
Jianhong Ou
}
\examples{
  path <- system.file("extdata", package = "motifStack")
  importMatrix(dir(path, "*.pcm", full.names = TRUE))
}
\keyword{ misc }
