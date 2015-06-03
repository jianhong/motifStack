\name{motifCircos}
\alias{motifCircos}
\title{
plot sequence logo stacks with a radial phylogenic tree and multiple color rings
}
\description{
plot sequence logo stacks with a radial phylogenic tree and multiple color rings.
The difference from plotMotifStackWithRadialPhylog is that it has more color setting
and one more group of pfms.
}
\usage{
motifCircos(phylog, pfms=NULL, pfms2=NULL, R=2.5,
            r.tree=1, col.tree.bg=NULL, col.tree.bg.alpha=1,
            cnodes=0, labels.nodes=names(phylog$nodes), clabel.nodes=0,
            r.leaves=NA,
            cleaves=1, labels.leaves=names(phylog$leaves), clabel.leaves=1,
            col.leaves=rep("black", length(labels.leaves)),
            col.leaves.bg=NULL, col.leaves.bg.alpha=1,
            r.pfms=NA, r.pfms2=NA,
            r.rings=0, col.rings=list(),
            col.inner.label.circle=NULL, inner.label.circle.width=0.02,
            col.outer.label.circle=NULL, outer.label.circle.width=0.02,
            draw.box=FALSE,
            clockwise =FALSE, init.angle=if(clockwise) 90 else 0,
            angle=360, pfmNameSpliter=";", rcpostfix="(RC)", 
            motifScale=c("linear","logarithmic"), ic.scale=TRUE,
            plotIndex=FALSE, IndexCol="black", IndexCex=.8,
            groupDistance=NA, groupDistanceLineCol="red", 
            plotAxis=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{phylog}{an object of class phylog}
\item{pfms}{a list of objects of class pfm}
\item{pfms2}{a list of objects of class pfm}
\item{R}{radius of canvas}
\item{r.tree}{half width of the tree}
\item{col.tree.bg}{a vector of colors for tree background}
\item{col.tree.bg.alpha}{a alpha value [0, 1] of colors for tree background}
\item{cnodes}{a character size for plotting the points that represent the nodes, used with par("cex")*cnodes. If zero, no points are drawn}
\item{labels.nodes}{a vector of strings of characters for the nodes labels}
\item{clabel.nodes}{a character size for the nodes labels, used with par("cex")*clabel.nodes. If zero, no nodes labels are drawn}
\item{r.leaves}{width of the leaves}
\item{cleaves}{a character size for plotting the points that represent the leaves, used with par("cex")*cleaves. If zero, no points are drawn}
\item{labels.leaves}{a vector of strings of characters for the leaves labels}
\item{clabel.leaves}{a character size for the leaves labels, used with}
\item{col.leaves}{a vector of colors for leaves labels}
\item{col.leaves.bg}{a vector of colors for background of leaves labels}
\item{col.leaves.bg.alpha}{alpha value [0, 1] for the colors of backgroud of leaves labels}
\item{r.pfms}{width of the pfms}
\item{r.pfms2}{width of the pfms2}
\item{r.rings}{a vector of width of color rings}
\item{col.rings}{a list of color rings}
\item{col.inner.label.circle}{a vector of colors for inner cirlce of pfms}
\item{inner.label.circle.width}{width for inner circle of pfms}
\item{col.outer.label.circle}{a vector of colors for outer circle of pfms}
\item{outer.label.circle.width}{width for outer circle of pfms}
\item{draw.box}{if TRUE draws a box around the current plot with the function box()}
\item{clockwise}{a logical value indicating if slices are drawn clockwise or counter clockwise}
\item{init.angle}{number specifying the starting angle (in degrees) for the slices. Defaults to 0 (i.e., `3 o'clock`) unless clockwise is true where init.angle defaults to 90 (degrees), (i.e., `12 o'clock`)}
\item{angle}{number specifying the angle (in degrees) for phylogenic tree. Defaults 360}
\item{pfmNameSpliter}{spliter when name of pfms/pfms2 contain multiple node of labels.leaves}
\item{rcpostfix}{the postfix for reverse complements}
\item{motifScale}{the scale of logo size}
\item{ic.scale}{logical. If TRUE, the height of each column is proportional to its information content. Otherwise, all columns have the same height.}
\item{plotIndex}{logical. If TRUE, will plot index number in the motifLogo which can help user to describe the motifLogo}
\item{IndexCol}{The color of the index number when plotIndex is TRUE.}
\item{IndexCex}{The cex of the index number when plotIndex is TRUE.}
\item{groupDistance}{show groupDistance on the draw}
\item{groupDistanceLineCol}{groupDistance line color, default: red}
\item{plotAxis}{logical. If TRUE, will plot distance axis.}
}
\value{
none
}
\author{
Jianhong Ou
}
\seealso{
\code{\link{plotMotifStackWithRadialPhylog}}
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
    motifCircos(phylog, pfms, cleaves = 0.5, clabel.leaves = 0.7, 
                     col.tree.bg=rep(color, each=5), col.leaves=rep(color, each=5),
                      r.rings=c(0.02, 0.03, 0.04), 
                      col.rings=list(sample(colors(), 50), 
                                     sample(colors(), 50), 
                                     sample(colors(), 50)))
  }
}
\keyword{ misc }
