\name{motifCloud}
\alias{motifCloud}
\title{plot a DNA sequence logo cloud}
\description{
  Plot a DNA sequence logo cloud
}
\usage{
motifCloud(motifSig, rcpostfix="(RC)", 
           layout=c("rectangles", "cloud", "tree"), 
           scale=c(6, .5), rot.per=.1,
           draw.box=TRUE, draw.freq=TRUE, 
           box.col="gray", freq.col="gray",
           group.col=NULL, groups=NULL, draw.legend=FALSE, 
           font="Helvetica-Bold", ic.scale=TRUE)
}
\arguments{
\item{motifSig}{an object of class \linkS4class{motifSig}}
\item{rcpostfix}{postfix for reverse-complement motif names, default: (RC)}
\item{layout}{layout of the logo cloud, rectangles, cloud or tree}
\item{scale}{A vector of length 2 indicating the range of the size of the sequence logo. }
\item{rot.per}{proportion sequence logo with 90 degree rotation. Only work for "cloud" layout}
\item{draw.box}{draw box for each sequence logo or not}
\item{draw.freq}{label frequency of each signature or not}
\item{box.col}{color of box for each sequence logo}
\item{freq.col}{color of frequency label}
\item{group.col}{color setting for groups}
\item{groups}{a named vectors of motif groups}
\item{draw.legend}{draw group color legend or not}
\item{font}{font of logo}
\item{ic.scale}{logical If TRUE, the height of each column is proportional to its information content. Otherwise, all columns have the same height.}
}
\value{
 none
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
    d <- MotIV::motifDistances(lapply(pfms, pfm2pwm))
    hc <- MotIV::motifHclust(d, method="average")
    phylog <- hclust2phylog(hc)
    leaves <- names(phylog$leaves)
    pfms <- pfms[leaves]
    pfms <- lapply(names(pfms), function(.ele, pfms){new("pfm",mat=pfms[[.ele]], 
                                                     name=.ele)},pfms)
    motifSig <- motifSignature(pfms, phylog, groupDistance=0.1)
    motifCloud(motifSig)
  }
}