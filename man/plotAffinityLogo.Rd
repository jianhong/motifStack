\name{plotAffinityLogo}
\alias{plotAffinityLogo}
\title{plot affinity logo}
\description{
plot amino acid or DNA sequence affinity logo
}
\usage{
plotAffinityLogo(psam, motifName, font="Helvetica-Bold",
                             colset=c("#00811B","#2000C7","#FFB32C","#D00001"),
                             alpha=0.5, newpage=TRUE, fontsize=12)
}
\arguments{
\item{psam}{a position frequency matrices}
\item{motifName}{motif name}
\item{font}{font of logo}
\item{colset}{color setting for each logo letter}
\item{alpha}{Alpha channel for transparency of low affinity letters.}
\item{newpage}{plot in a new canvas or not.}
\item{fontsize}{font size of the template for grImport, default 12. Higher value make better quality figure, but also increase the file size.}
}
\value{
 none
}
\examples{
psam <- importMatrix(file.path(find.package("motifStack"), "extdata", "PSAM.mxr"), 
                     format="psam")[[1]]
plotAffinityLogo(psam)
}