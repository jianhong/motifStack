\name{plotMotifStackWithRadialPhylog}
\alias{plotMotifStackWithRadialPhylog}
\title{plot sequence logo stacks with a radial phylogenic tree}
\description{
plot sequence logo stacks with a radial phylogenic tree
}
\usage{
plotMotifStackWithRadialPhylog(phylog, pfms=NULL,  
circle=1, cleaves=1, cnodes=0, 
labels.leaves=names(phylog$leaves), clabel.leaves=1,
labels.nodes=names(phylog$nodes), clabel.nodes=0, 
draw.box=FALSE, 
col.leaves=rep("black", length(labels.leaves)), 
col.leaves.bg=NULL, col.leaves.bg.alpha=1,
col.bg=NULL, col.bg.alpha=1, 
col.inner.label.circle=NULL, col.outer.label.circle=NULL, outer.label.circle.width="default",
clockwise =FALSE, init.angle=if(clockwise) 90 else 0, 
angle=360)
}
\arguments{
\item{phylog}{an object of class phylog}
\item{pfms}{a list of objects of class pfm}
\item{circle}{a size coefficient for the outer circle}
\item{cleaves}{a character size for plotting the points that represent the leaves, used with par("cex")*cleaves. If zero, no points are drawn}
\item{cnodes}{a character size for plotting the points that represent the nodes, used with par("cex")*cnodes. If zero, no points are drawn}
\item{labels.leaves}{a vector of strings of characters for the leaves labels}
\item{clabel.leaves}{a character size for the leaves labels, used with}
\item{labels.nodes}{a vector of strings of characters for the nodes labels}
\item{clabel.nodes}{a character size for the nodes labels, used with par("cex")*clabel.nodes. If zero, no nodes labels are drawn}
\item{draw.box}{if TRUE draws a box around the current plot with the function box()}
\item{col.leaves}{a vector of colors for leaves labels}
\item{col.leaves.bg}{a vector of colors for background of leaves labels}
\item{col.leaves.bg.alpha}{alpha value [0, 1] for the colors of backgroud of leaves labels}
\item{col.bg}{a vector of colors for background}
\item{col.bg.alpha}{a alpha value [0, 1] of colors for background}
\item{col.inner.label.circle}{a vector of colors for inner cirlce of leaves labels}
\item{col.outer.label.circle}{a vector of colors for outer circle of leaves labels}
\item{outer.label.circle.width}{width for outer circle of leaves labels}
\item{clockwise}{a logical value indicating if slices are drawn clockwise or counter clockwise}
\item{init.angle}{number specifying the starting angle (in degrees) for the slices. Defaults to 0 (i.e., `3 o'clock`) unless clockwise is true where init.angle defaults to 90 (degrees), (i.e., `12 o'clock`)}
\item{angle}{number specifying the angle (in degrees) for phylogenic tree. Defaults 360}
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
		library("MotIV")
		matrix.fly <- query(MotifDb, "Dmelanogaster")
		motifs <- as.list(matrix.fly)
		motifs <- motifs[grepl("Dmelanogaster-FlyFactorSurvey-", names(motifs), fixed=TRUE)]
		names(motifs) <- gsub("Dmelanogaster_FlyFactorSurvey_", "", 
							  gsub("_FBgn[0-9]+$", "", 
									gsub("[^a-zA-Z0-9]","_", 
										 gsub("(_[0-9]+)+$", "", names(motifs)))))
		motifs <- motifs[unique(names(motifs))]
		pfms <- sample(motifs, 50)
		jaspar <- MotIV::readPWMfile(file.path(find.package("MotIV"), "extdata", "jaspar2010.txt"))
		jaspar.scores <- MotIV::readDBScores(file.path(find.package("MotIV"), "extdata", "jaspar2010_PCC_SWU.scores"))
		ffs.analysis.jaspar<-MotIV::motifMatch(inputPWM=pfms,align="SWU",cc="PCC",database=jaspar, DBscores=jaspar.scores,top=5)
		d <- MotIV::motifDistances(getPWM(ffs.analysis.jaspar))
		hc <- MotIV::motifHclust(d)
		phylog <- hclust2phylog(hc)
		leaves <- names(phylog$leaves)
		pfms <- pfms[leaves]
		pfms <- lapply(names(pfms), function(.ele, pfms){new("pfm",mat=pfms[[.ele]], name=.ele)},pfms)
		pfms <- DNAmotifAlignment(pfms, minimalConsensus=3)
		library(RColorBrewer)
		color <- brewer.pal(12, "Set3")
		plotMotifStackWithRadialPhylog(phylog, pfms, circle=0.9, cleaves = 0.5, clabel.leaves = 0.7, 
										 col.bg=rep(color, each=5), col.leaves=rep(color, each=5))
	}
}