#' Matrix Aligner
#' @description Matrix Aligner is modified from Matalign-v4a.
#' Matalign-v4a is a program to compare two positional specific matrices.
#' The author of Matalign-v4a is Ting Wang and Gary Stormo.
#' @param pcms A list of \link{pcm}
#' @param method Alignment method. "Smith-Waterman" or "Needleman-Wunsch".
#' Default is "Smith-Waterma"
#' @param pseudo pseudocount
#' @param revComp Check reverseComplement or not.
#' @return A data frame with alignment information. The column names
#' are motif1, motif2, alignmentScore, startPos1, startPos2, endPos1, endPos2,
#' alignmentLength.
#' @export
#' @importFrom utils combn
#' @examples 
#' if(interactive() || Sys.getenv("USER")=="jianhongou"){
#'   fp <- system.file("extdata", package="motifStack")
#'   fs <- dir(fp, "pcm$")
#'   pcms <- importMatrix(file.path(fp, fs), format="pcm")
#'   matalign(pcms)
#'  }
matalign <- function(pcms, 
                     method=c("Smith-Waterman", "Needleman-Wunsch"),
                     pseudo=1,
                     revComp=TRUE){
  method <- match.arg(method, 
                      choices = c("Smith-Waterman", "Needleman-Wunsch"))
  stopifnot(is.numeric(pseudo))
  stopifnot(is.list(pcms))
  if(length(names(pcms))!=length(pcms)) 
    names(pcms) <- paste0("M", seq_along(pcms))
  pcms <- mapply(pcms, names(pcms), FUN=function(.ele, .name){
    if(is(.ele, "pcm")){
      return(.ele)
    }
    if(is(.ele, "pfm")){
      return(pfm2pcm(.ele))
    }
    if(is(.ele, "matrix") || is(.ele, "data.frame")){
      cols <- colSums(.ele)
      if(all(abs(1 - cols) < 0.01)){
        .ele <- new("pfm", name=.name, mat=.ele)
        return(pfm2pcm(.ele))
      }
      if(all(cols>=1) && all(cols==round(cols))){
        return(new("pcm", name=.name, mat=.ele))
      }
    }
    stop("motifs must be pcms")
  }, SIMPLIFY = FALSE)
  stopifnot("pcms must be a list of pcm object" = is.list(pcms))
  null <- lapply(pcms, function(.ele) if(!is(.ele, "pcm")){
    stop("pcms must be a list of pcm object.")
  })
  n <- vapply(pcms, FUN = function(.ele) .ele@name,
              FUN.VALUE = "a", USE.NAMES = FALSE)
  stopifnot("There is duplicated motif names"=!duplicated(n))
  names(pcms) <- n
  cmb <- combn(n, m = 2, simplify = FALSE)
  align <- lapply(cmb, FUN=function(.ele){
    compareProfiles(pcms[[.ele[1]]], pcms[[.ele[2]]],
                    pseudo=pseudo, revComp=revComp)
  })
  align <- do.call(rbind, align)
  align <- cbind(do.call(rbind, cmb), align[, -(2:3)])
  colnames(align) <- c("motif1", "motif2", "alignmentScore", 
                      "startPos1", "startPos2", "endPos1", "endPos2",
                      "alignmentLength", "P_value", 
                      "distance", 'alignedDist',
                      "direction")[seq.int(ncol(align))]
  return(align)
}

#' comapre w pcm
#' @description compare two pcm objects
#' @param pcm1,pcm2 object of pcm
#' @param method Alignment method. "Smith-Waterman" or "Needleman-Wunsch".
#' Default is "Smith-Waterma"
#' @param pseudo pseudocount
#' @param revComp Check reverseComplement or not.
#' @return a list with names: motif1, motif2, alignmentScore, 
#' startPos1, startPos2, endPos1, endPos2, alignmentLength.
compareProfiles <- function(pcm1, pcm2, 
                            method=c("Smith-Waterman", "Needleman-Wunsch"),
                            pseudo=1, 
                            revComp=TRUE){
  method <- match.arg(method, 
                      choices = c("Smith-Waterman", "Needleman-Wunsch"))
  hspF <- compare2profiles(pcm1, pcm2, method, pseudo)
  if(!(pcm1@alphabet %in% c("DNA", "RNA")) || !revComp){
    hspF$direction <- "Forward"
    return(hspF)
  }
  hspR <- compare2profiles(pcm1, matrixReverseComplement(pcm2),
                           method, pseudo)
  if(hspF$score >= hspR$score){
    hspF$direction <- "Forward"
    return(hspF)
  }else{
    hspR$direction <- "Reverse"
    return(hspR)
  }
}

#' compare two profiles
#' @description compare two pcm object
#' @param pcm1,pcm2 object of pcm
#' @param method Alignment method. "Smith-Waterman" or "Needleman-Wunsch".
#' Default is "Smith-Waterma"
#' @param pseudo pseudocount
#' @return a list with names: motif1, motif2, alignmentScore, 
#' startPos1, startPos2, endPos1, endPos2, alignmentLength.
compare2profiles <- function(pcm1, pcm2, 
                             method = c("Smith-Waterman", 
                                        "Needleman-Wunsch"),
                             pseudo=1
                             ){
  method <- match.arg(method, 
                      choices = c("Smith-Waterman", "Needleman-Wunsch"))
  profile1 <- pcm1@mat
  profile2 <- pcm2@mat
  width1 <- ncol(profile1)
  width2 <- ncol(profile2)
  
  score <- getScore(pcm1, pcm2, pseudo)
  
  if(method == "Needleman-Wunsch"){
    dpScore <- dpGlobal(score, width1, width2)
    hsp <- traceBackGlobal(dpScore, score, width1, width2)
  }else{
    dpScore <- dpLocal(score, width1, width2)
    hsp <- traceBackLocal(dpScore, score, width1, width2)
  }
  dist <- getDistance(hsp, count1 = profile1, count2 = profile2,
                      P1 = pcm1@background, P2 = pcm2@background,
                      pseudo = pseudo)
  hsp$distance <- dist[1]
  hsp$alignedDist <- dist[2]
  return(hsp)
}


#' Calculate pair_wise position score
#' @param pcm1,pcm2 object of pcm
#' @param pseudo pseudocount
#' @return A score matrix with nrow of ncol of pcm1 
#' and ncol of ncol of pcm2.
getScore <- function(pcm1, pcm2, pseudo=1){
  stopifnot("input is not pcm object"=is(pcm1, "pcm"))
  stopifnot("input is not pcm object"=is(pcm2, "pcm"))
  profile1 <- pcm1@mat
  profile2 <- pcm2@mat
  P1 <- pcm1@background
  P2 <- pcm2@background
  score <- matrix(0, 
                  nrow = ncol(profile1),
                  ncol = ncol(profile2))
  for(i in seq.int(ncol(profile1))){
    for(j in seq.int(ncol(profile2))){
      score[i, j] <- 
        getALLRscoreFromCounts(profile1[, i], profile2[, j], 
                               P1, P2, pseudo)
    }
  }
  return(score)
}

#' calculate I'
#' @param freq1 position frequence for matrix 1 position j
#' @param freq2 position frequence for matrix 2 position j
#' @param P background of profile1
#' @return numeric(1)
calI <- function(freq1, freq2, P){
  sum(freq2*log(freq1/P)/log(2))
}

#' calculate frequence
#' @param count position counts
#' @param P background probility 
#' @param pseudo pseudocount
#' @return numeric(1)
calF <- function(count, 
                 P=rep(1/length(count), length(count)),
                 pseudo=1){
  total <- sum(count)
  f <- (count + P*pseudo)/(total+pseudo)
  return(f)
}

#' calculate ALLR from counts
#' @param count1,count2 count in position j for matrix 1 or 2
#' @param P1,P2 background for matrix 1 or 2
#' @param pseudo pseudocount
#' @return numeric(1) of ALLR
getALLRscoreFromCounts <- function(count1, count2, P1, P2, pseudo){
  freq1 <- calF(count1, P1, pseudo)
  freq2 <- calF(count2, P2, pseudo)
  
  I1 <- calI(freq1, freq2, P1)
  I2 <- calI(freq2, freq1, P2)
  
  total1 <- sum(count1)
  total2 <- sum(count2)
  
  score <- (I1*total2 + I2*total1)/(total1+total2)
  return(score)
}

#' Dynamic programming function, local version 
#' @param score ALLR scores, m x n matrix
#' @param m,n matrix width 
#' @return score matrix
dpLocal <- function(score, m, n){
  s <- 0
  DP <- matrix(0, nrow = m, ncol = n)
  DP[score[, 1]>0, 1] <- score[score[, 1]>0, 1]
  DP[1, score[1, ]>0] <- score[1, score[1, ]>0]
  for(i in seq.int(m)[-1]){
    for(j in seq.int(n)[-1]){
      s <- DP[i-1, j-1] + score[i, j]
      DP[i, j] <- ifelse(s>0, s, 0)
    }
  }
  return(DP)
}

#' Global alignment version
#' @param score ALLR scores, m x n matrix
#' @param m,n matrix width 
#' @return score matrix
dpGlobal <- function(score, m, n){
  s <- 0
  DP <- matrix(0, nrow = m, ncol = n)
  DP[, 1] <- score[, 1]
  DP[1, ] <- score[1, ]
  for(i in seq.int(m)[-1]){
    for(j in seq.int(n)[-1]){
      s <- DP[i-1, j-1] + score[i, j]
      DP[i, j] <- s
    }
  }
  return(DP)
}

#' traceback local
#' @param dpScore Dynamic programming score matrix
#' @param score ALLR scores, m x n matrix
#' @param m,n matrix width
#' @return a data.frame
traceBackLocal <- function(dpScore, score, m, n){
  max <- max(dpScore)
  maxpos <- which.max(dpScore)
  maxpos <- maxpos[length(maxpos)]
  j <- j_end <- ceiling(maxpos/m)
  i <- i_end <- maxpos - m*(j_end - 1)
  while(i>1 && j>1){
    if(dpScore[i, j]>0){
      i <- i-1
      j <- j-1
    }else{
      i <- i+1
      j <- j+1
      break()
    }
  }
  i_start <- i
  j_start <- j
  eval <- getEval(max, m, n)
  pval <- Eval2Pval(eval)
  return(data.frame(score=max,
              width1=m,
              width2=n,
              i_start=i_start,
              j_start=j_start,
              i_end=i_end,
              j_end=j_end,
              length=i_end - i_start + 1,
              P_value=pval))
}
#' traceback global
#' @param dpScore Dynamic programming score
#' @param score ALLR scores
#' @param m,n matrix width 
#' @return a data.frame
traceBackGlobal <- function(dpScore, score, m, n){
  max <- dpScore[m, n]
  i_end <- m
  j_end <- n
  for(i in seq.int(m)){
    if(dpScore[i, n]>=max){
      max <- dpScore[i, n]
      i_end <- i
      j_end <- n
    }
  }
  for(j in seq.int(n)){
    if(dpScore[m, j]>=max){
      max <- dpScore[m, j]
      i_end <- m
      j_end <- j
    }
  }
  i <- i_end
  j <- j_end
  while(i>1 && j>1){
    i <- i-1
    j <- j-1
  }
  i_start <- i
  j_start <- j
  eval <- getEval(max, m, n)
  pval <- Eval2Pval(eval)
  return(data.frame(score=max,
                    width1=m,
                    width2=n,
                    i_start=i_start,
                    j_start=j_start,
                    i_end=i_end,
                    j_end=j_end,
                    length=i_end - i_start + 1,
                    P_value=pval))
}

#' Calculate distances between two profiles
#' @param hsp output of traceBack function
#' @param count1,count2 motif profile 1 or 2
#' @param P1,P2 background of profile 1 or 2
#' @param pseudo pseudocount
#' @return full distance and aligned distance.
getDistance <- function(hsp, count1, count2, 
                        P1, P2, pseudo){
  full1 <- sum(
    apply(count1, 2, function(.ele) 
      getALLRscoreFromCounts(count1 = .ele, count2 = .ele, 
                             P1=P1, P2=P1, pseudo = pseudo))
  )
  if(hsp$i_start<=hsp$i_end){
    aligned1 <- sum(
      apply(count1[, hsp$i_start:hsp$i_end, drop=FALSE], 2, function(.ele) 
        getALLRscoreFromCounts(count1 = .ele, count2 = .ele, 
                               P1=P1, P2=P1, pseudo = pseudo))
    )
  }else{
    aligned1 <- 0
  }
  
  full2 <- sum(
    apply(count2, 2, function(.ele) 
      getALLRscoreFromCounts(count1 = .ele, count2 = .ele, 
                             P1=P2, P2=P2, pseudo = pseudo))
  )
  if(hsp$i_start<=hsp$i_end){
    aligned2 <- sum(
      apply(count2[, hsp$j_start:hsp$j_end, drop=FALSE], 2, function(.ele) 
        getALLRscoreFromCounts(count1 = .ele, count2 = .ele,
                               P1=P2, P2=P2, pseudo = pseudo))
    )
  }else{
    aligned2 <- 0
  }
  return(c(
    dist1 = full1 + full2 - 2*hsp$score,
    dist2 = aligned1 + aligned2 - 2*hsp$score
  ))
}

getEval <- function(alignmentScore, m, n){
  Lambda <- 0.02094
  K <- 0.26062
  H <- 3.25810
  K * m * n * exp (-Lambda*(alignmentScore*100))
}

Eval2Pval <- function(e){
  1 - exp(-e)
}
