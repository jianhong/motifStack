\name{plotMotifStackWithRadialPhylog}
\alias{plotMotifStackWithRadialPhylog}
\title{plot sequence logo stacks with a radial phylogenic tree}
\description{
plot sequence logo stacks with a radial phylogenic tree
}
\usage{
plotMotifStackWithRadialPhylog(phylog, pfms=NULL,  
circle = 1, cleaves = 1, cnodes = 0, 
labels.leaves = names(phylog$leaves), clabel.leaves = 1, 
labels.nodes = names(phylog$nodes), clabel.nodes = 0, 
draw.box = FALSE, col.leaves=rep("black", length(labels.leaves)),
col.bg=NULL)
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
\item{col.bg}{a vector of colors for background}
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
		motifs <- motifs[grepl("Dmelanogaster\\-FlyFactorSurvey\\-", names(motifs))]
		names(motifs) <- gsub("Dmelanogaster_FlyFactorSurvey_", "", 
							  gsub("_FBgn\\d+$", "", 
									gsub("[^a-zA-Z0-9]","_", 
										 gsub("(_\\d+)+$", "", names(motifs)))))
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