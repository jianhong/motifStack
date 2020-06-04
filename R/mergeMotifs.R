#' merge multiple motifs
#' 
#' merge multiple motifs by calculate mean of each position
#' 
#' 
#' @param ... \link{pcm} or \link{pfm} objects
#' @param bgNoise if it is not NA, test will using a background by
#' Dirichlet(1)-distributed random frequencies with weight bg.noise.  The value
#' of bgNoise should be a number in the range of 0 to 1, eg. 0.05
#' @return a \link{pfm} object
#' @author Jianhong Ou
#' @keywords misc
#' @export
#' @importFrom stats rgamma
#' @examples
#' 
#'     pcms<-readPCM(file.path(find.package("motifStack"), "extdata"),"pcm$")
#'     mergeMotifs(pcms)
#' 
mergeMotifs <- function(..., bgNoise=NA){
    mL <- list(...)
    n <- length(mL)
    if(n==1){
        mL <- mL[[1]]
        n <- length(mL)
    }
    if(n<2){
        stop("The number of input motifs should more than 1")
    }
    classes <- sapply(mL, function(.ele) {
      if(is(.ele, "pcm")) return("pcm")
      if(is(.ele, "pfm")) return("pfm")
      return("others")
    })
    if(any(!classes %in% c("pcm", "pfm"))){
        stop("Not all of the motifs are pcm or pfm object")
    }
    if(length(unique(classes))!=1){
        stop("The inputs are not in same format")
    }
    names <- sapply(mL, function(.ele) .ele$name)
    if(length(names)!=n){
        message("Not all motifs have names. the Name slot will be set as NULL")
        names <- ""
    }else{
        names <- paste(names, collapse=";")
    }
    
    mL <- DNAmotifAlignment(mL)
    d <- sapply(mL, function(.ele) dim(.ele$mat))
    if(nrow(unique(t(d)))!=1) stop("Unexpect happened. Report this bug please.")
    d <- c(d[, 1], n)
    dat <- lapply(mL, function(.ele) .ele$mat)
    dat <- array(unlist(dat), dim=d)
    mat <- apply(dat, MARGIN=c(1, 2), FUN=mean)
    mat <- t(t(mat)/colSums(mat)) #convert into frequency
    
    rdirichlet <- function (n, alpha){
        l <- length(alpha)
        x <- matrix(rgamma(l * n, alpha), ncol = n, byrow = TRUE)
        sm <- x %*% rep(1, n)
        t(x/as.vector(sm))
    }
    if(!is.na(bgNoise)){
        if(bgNoise<0 | bgNoise>1) stop("bgNoise must be number in (0, 1)")
        mat <- (1-bgNoise)*mat + 
            bgNoise*rdirichlet(nrow(mat), rep(1, ncol(mat)))
    }
    
    rownames(mat) <- rownames(mL[[1]]$mat)
    
    new("pfm", mat=mat, name=names)
}
