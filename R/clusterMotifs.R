clusterMotifs <- function(pfms){
  if(requireNamespace("MotIV", quietly = TRUE)){
    jaspar.scores <- 
      MotIV::readDBScores(file.path(find.package("MotIV"), 
                                    "extdata", 
                                    "jaspar2010_PCC_SWU.scores"))
    d <- MotIV::motifDistances(pfmList2matrixList(pfms), 
                               DBscores=jaspar.scores)
    hc <- MotIV::motifHclust(d, method="average")
  }else{
    stop("MotIV is required.")
  }
  hc
}


pfmList2matrixList <- function(pfms){
  m <- lapply(pfms, pfm2pwm)
  names(m) <- unlist(lapply(pfms, function(.ele) .ele@name))
  m
}