\name{motifCloud}
\alias{motifCloud}
\title{plot a DNA sequence logo cloud}
\description{
	Plot a DNA sequence logo cloud
}
\usage{
motifCloud(pfms, phylog, groupDistance, rcpostfix="(RC)", 
layout=c("rectangles", "cloud", "tree"), 
min.freq=2, scale=c(6, .5), rot.per=.1, trim=0.8, 
draw.box=TRUE, draw.freq=TRUE, 
box.col="gray", freq.col="gray",
group.col=NULL, groups=NULL, draw.legend=FALSE)
}
\arguments{
\item{pfms}{a list of objects of class pfm}
\item{phylog}{an object of class phylog}
\item{groupDistance}{maxmal distance of motifs in the same group}
\item{rcpostfix}{postfix for reverse-complement motif names, default: (RC)}
\item{layout}{layout of the logo cloud, rectangles, cloud or tree}
\item{min.freq}{signatures with frequency below min.freq will not be plotted}
\item{scale}{A vector of length 2 indicating the range of the size of the sequence logo. }
\item{rot.per}{proportion sequence logo with 90 degree rotation. Only work for "cloud" layout}
\item{trim}{minimal information content for each position of signature}
\item{draw.box}{draw box for each sequence logo or not}
\item{draw.freq}{label frequency of each signature or not}
\item{box.col}{color of box for each sequence logo}
\item{freq.col}{color of frequency label}
\item{group.col}{color setting for groups}
\item{groups}{a named vectors of motif groups}
\item{draw.legend}{draw group color legend or not}
}
\value{
 none
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
		motifCloud(pfms, phylog, groupDistance=0.1)
	}
}