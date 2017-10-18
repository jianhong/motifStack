\name{motifPiles}
\alias{motifPiles}
\title{
plot sequence logo stacks with a linear phylogenic tree and multiple color sets
}
\description{
plot sequence logo stacks with a linear phylogenic tree and multiple color sets.
}
\usage{
motifPiles(phylog, pfms=NULL, pfms2=NULL, 
           r.tree=.45, col.tree=NULL,
           cnodes=0, labels.nodes=names(phylog$nodes), clabel.nodes=0,
           cleaves=.2, labels.leaves=names(phylog$leaves), clabel.leaves=1,
           col.leaves=rep("black", length(labels.leaves)),
           col.leaves.bg=NULL, col.leaves.bg.alpha=1,
           r.pfms=NA, r.pfms2=NA, motifScale=c("logarithmic", "linear", "none"),
           col.pfms=NULL, col.pfms.width=0.02,
           col.pfms2=NULL, col.pfms2.width=0.02,
           r.anno=0, col.anno=list(),
           pfmNameSpliter=";", rcpostfix="(RC)", ic.scale=TRUE,
           plotIndex=FALSE, IndexCol="black", IndexCex=.8,
           groupDistance=NA, groupDistanceLineCol="red")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{phylog}{an object of class phylog}
\item{pfms}{a list of objects of class pfm}
\item{pfms2}{a list of objects of class pfm}
\item{r.tree}{width of the tree}
\item{col.tree}{a vector of colors for tree}
\item{cnodes}{a character size for plotting the points that represent the nodes, used with par("cex")*cnodes. If zero, no points are drawn}
\item{labels.nodes}{a vector of strings of characters for the nodes labels}
\item{clabel.nodes}{a character size for the nodes labels, used with par("cex")*clabel.nodes. If zero, no nodes labels are drawn}
\item{cleaves}{a character size for plotting the points that represent the leaves, used with par("cex")*cleaves. If zero, no points are drawn}
\item{labels.leaves}{a vector of strings of characters for the leaves labels}
\item{clabel.leaves}{a character size for the leaves labels, used with}
\item{col.leaves}{a vector of colors for leaves labels}
\item{col.leaves.bg}{a vector of colors for background of leaves labels}
\item{col.leaves.bg.alpha}{alpha value [0, 1] for the colors of backgroud of leaves labels}
\item{r.pfms}{width of the pfms}
\item{r.pfms2}{width of the pfms2}
\item{motifScale}{the scale of logo size}
\item{col.pfms}{a vector of colors for inner pile of pfms}
\item{col.pfms.width}{width for inner pile of pfms}
\item{col.pfms2}{a vector of colors for outer pile of pfms}
\item{col.pfms2.width}{width for outer pile of pfms}
\item{r.anno}{a vector of width of color sets}
\item{col.anno}{a list of color sets}
\item{pfmNameSpliter}{spliter when name of pfms/pfms2 contain multiple node of labels.leaves}
\item{rcpostfix}{the postfix for reverse complements}
\item{ic.scale}{logical. If TRUE, the height of each column is proportional to its information content. Otherwise, all columns have the same height.}
\item{plotIndex}{logical. If TRUE, will plot index number in the motifLogo which can help user to describe the motifLogo}
\item{IndexCol}{The color of the index number when plotIndex is TRUE.}
\item{IndexCex}{The cex of the index number when plotIndex is TRUE.}
\item{groupDistance}{show groupDistance on the draw}
\item{groupDistanceLineCol}{groupDistance line color, default: red}
}
\value{
none
}
\author{
Jianhong Ou
}
\seealso{
\code{\link{motifCircos}}
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
    library(RColorBrewer)
    color <- brewer.pal(12, "Set3")
    motifPiles(phylog, pfms, cleaves = 0.5, clabel.leaves = 0.7, 
             col.leaves=rep(color, each=5), 
             col.leaves.bg = sample(colors(), 50),
             col.tree=rep(color, each=5),
             r.anno=c(0.02, 0.03, 0.04), 
             col.anno=list(sample(colors(), 50), 
                            sample(colors(), 50), 
                            sample(colors(), 50)))
  }
}
\keyword{ misc }
