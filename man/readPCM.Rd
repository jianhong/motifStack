\name{readPCM}
\alias{readPCM}
\title{read pcm from a path}
\description{
	read position count matrix from a path
}
\usage{
readPCM(path=".", pattern=NULL)
}
\arguments{
\item{path}{a character vector of full path names}
\item{pattern}{an optional regular expression}
}
\value{
	A list of \code{\link{pcm}} objects
}
\examples{
pcms<-readPCM(file.path(find.package("motifStack"), "extdata"),"pcm$")
}