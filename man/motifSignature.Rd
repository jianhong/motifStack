% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motifSignature.R
\name{motifSignature}
\alias{motifSignature}
\title{get signatures from motifs}
\usage{
motifSignature(
  pfms,
  phylog,
  cutoffPval,
  groupDistance,
  rcpostfix = "(RC)",
  min.freq = 2,
  trim = 0.2,
  families = list(),
  sort = TRUE
)
}
\arguments{
\item{pfms}{a list of objects of class pfm}

\item{phylog}{an object of class phylog}

\item{cutoffPval}{pvalue for motifs to merge.}

\item{groupDistance}{maxmal distance of motifs in the same group}

\item{rcpostfix}{postfix for reverse-complement motif names, default: (RC)}

\item{min.freq}{signatures with frequency below min.freq will not be plotted}

\item{trim}{minimal information content for each position of signature}

\item{families}{for each family, the motif number in one signature should
only count as 1}

\item{sort}{sort the signatures by frequency or not.}
}
\value{
an Object of class \linkS4class{motifSig}
}
\description{
extract signatures from multiple motifs by distance calculated from STAMP
}
\examples{

  if(interactive() || Sys.getenv("USER")=="jianhongou"){
    library("MotifDb")
    matrix.fly <- query(MotifDb, "Dmelanogaster")
    motifs <- as.list(matrix.fly)
    motifs <- motifs[grepl("Dmelanogaster-FlyFactorSurvey-", 
                            names(motifs), fixed=TRUE)]
    names(motifs) <- gsub("Dmelanogaster_FlyFactorSurvey_", "", 
                gsub("_FBgn[0-9]+$", "", 
                  gsub("[^a-zA-Z0-9]","_", 
                     gsub("(_[0-9]+)+$", "", names(motifs)))))
    motifs <- motifs[unique(names(motifs))]
    pfms <- sample(motifs, 50)
    hc <- clusterMotifs(pfms)
    library(ade4)
    phylog <- ade4::hclust2phylog(hc)
    leaves <- names(phylog$leaves)
    pfms <- pfms[leaves]
    pfms <- mapply(pfms, names(pfms), FUN=function(.ele, .name){
                 new("pfm",mat=.ele, name=.name)})
    motifSig <- motifSignature(pfms, phylog, cutoffPval=0.0001)
  }

}
