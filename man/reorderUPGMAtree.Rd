\name{reorderUPGMAtree}
\alias{reorderUPGMAtree}
\title{
    re-order UPGMA tree
}
\description{
    re-order the UPGMA tree by adjacent motif distance
}
\usage{
reorderUPGMAtree(phylog, motifs, rcpostfix = "(RC)")
}
\arguments{
  \item{phylog}{an object of phylog}
  \item{motifs}{a list of objects of pfm}
  \item{rcpostfix}{the postfix for reverse complements}
}
\value{
    an object of phylog
}
\author{
    Jianhong Ou
}
\examples{
  if(interactive()){
    library("MotifDb")
    matrix.fly <- query(MotifDb, "Dmelanogaster")
    motifs <- as.list(matrix.fly)
    motifs <- motifs[grepl("Dmelanogaster-FlyFactorSurvey-", names(motifs), fixed=TRUE)]
    names(motifs) <- gsub("Dmelanogaster_FlyFactorSurvey_", "", 
                gsub("_FBgn[0-9]+$", "", 
                  gsub("[^a-zA-Z0-9]","_", 
                     gsub("(_[0-9]+)+$", "", names(motifs)))))
    motifs <- motifs[unique(names(motifs))]
    pfms <- sample(motifs, 50)
    jaspar.scores <- MotIV::readDBScores(file.path(find.package("MotIV"), 
                                   "extdata", "jaspar2010_PCC_SWU.scores"))
    d <- MotIV::motifDistances(pfms)
    hc <- MotIV::motifHclust(d, method="average")
    phylog <- hclust2phylog(hc)
    reorderUPGMAtree(phylog, pfms)
  }
}
\keyword{ misc }
