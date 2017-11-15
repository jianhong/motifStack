\name{plotMotifStackWithPhylog}
\alias{plotMotifStackWithPhylog}
\title{plot sequence logo stacks with a ape4-style phylogenic tree}
\description{
plot sequence logo stacks with a ape4-style phylogenic tree
}
\usage{
plotMotifStackWithPhylog(phylog, pfms=NULL,
                         f.phylog = 0.3, f.logo = NULL, cleaves =1, cnodes =0,
                         labels.leaves = names(phylog$leaves), clabel.leaves=1,
                         labels.nodes = names(phylog$nodes), clabel.nodes = 0, 
                         font="Helvetica-Bold", ic.scale=TRUE, fontsize=12)
}
\arguments{
\item{phylog}{an object of class phylog}
\item{pfms}{a list of objects of class pfm}
\item{f.phylog}{a size coefficient for tree size (a parameter to draw the tree in proportion to leaves label)}
\item{f.logo}{a size coefficient for the motif}
\item{cleaves}{a character size for plotting the points that represent the leaves, used with par("cex")*cleaves. If zero, no points are drawn}
\item{cnodes}{a character size for plotting the points that represent the nodes, used with par("cex")*cnodes. If zero, no points are drawn}
\item{labels.leaves}{a vector of strings of characters for the leaves labels}
\item{clabel.leaves}{a character size for the leaves labels, used with par("cex")*clavel.leaves}
\item{labels.nodes}{a vector of strings of characters for the nodes labels}
\item{clabel.nodes}{a character size for the nodes labels, used with par("cex")*clabel.nodes. If zero, no nodes labels are drawn}
\item{font}{font of logo}
\item{ic.scale}{logical If TRUE, the height of each column is proportional to its information content. Otherwise, all columns have the same height.}
\item{fontsize}{font size of the template for grImport, default 12. Higher value make better quality figure, but also increase the file size.}
}
\value{
 none
}
\seealso{
\link[ade4:plot.phylog]{plot.phylog}
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
    pfms <- DNAmotifAlignment(pfms, minimalConsensus=3)
    plotMotifStackWithPhylog(phylog, pfms, f.phylog=0.3, 
                             cleaves = 0.5, clabel.leaves = 0.7)
  }
}