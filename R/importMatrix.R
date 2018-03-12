importMatrix <- function(filenames, 
                         format=c("auto", "pfm", "cm", "pcm", "meme", 
                                  "transfac", "jaspar", "scpd", "cisbp",
                                  "psam"), 
                         to=c("auto", "pcm", "pfm", "psam")){
  if(missing(filenames)){
    stop("filenames are required.")
  }
  stopifnot(is.character(filenames))
  format <- match.arg(format)
  to <- match.arg(to)
  stopifnot(all(file.exists(filenames)))
  fi <- file.info(filenames)
  if(any(fi$isdir)){
    stop("filenames could not contain folders.")
  }
  if(format=="auto"){
     getFileNameExtension <- function(fn){
       fn <- basename(fn)
       fn <- strsplit(fn, "\\.")
       ext <- lapply(fn, function(.ele){
         l <- length(.ele)
         if(l==1){
           NA
         }else{
           .ele[l]
         }
       })
       unlist(ext)
     }
     ext <- getFileNameExtension(filenames)
     if(any(is.na(ext))){
       stop("Can not determine the format of inputs by its extensions.")
     }
     ext <- unique(ext)
     if(length(ext)!=1){
       stop("There are multiple file formats in your inputs.")
     }
     format <- c("pfm", "cm", "pcm", "meme", "transfac", "jaspar", "scpd", "cisbp", "psam")
     format <- format[format==ext]
     if(length(format)!=1){
       stop("Can not determine the format of inputs by its extensions.")
     }
  }
  ## importM_... will output a list of matrix, 
  ## matrix rownames must be symbols, such as c(A, C, G, T)
  m <- do.call(paste0("importM_", format), list(fns=filenames))
  if(format=="meme"){
    return(m)
  }
  if(format=="psam"){
    return(m)
  }
  ## check is counts
  isCounts <- function(cnt){
    if(typeof(cnt)=="integer"){
      return(TRUE)
    }else{
      cnt1 <- round(cnt, digits = 0)
      if(all(cnt==cnt1)){
        return(TRUE)
      }
    }
    return(FALSE)
  }
  isCnt <- sapply(m, isCounts)
  # check if position frequency matrix
  isPFMmatrix <- function(cnt){
    if (any(abs(1 - colSums(cnt)) > 0.01)) return(FALSE)
    return(TRUE)
  }
  if(to=="auto"){
    if(all(isCnt)){
      to <- "pcm"
    }else{
      to <- "pfm"
    }
  }
  mapply(function(.m, .type, .f){
    if(to=="pcm"){
      if(!.type){
        stop("Can not import", .f, "to pcm format!")
      }
      new("pcm", mat=as.matrix(.m), name=.f)
    }else{
      if(isPFMmatrix(.m)){
        new("pfm", mat=as.matrix(.m), name=.f)
      }else{
        if(.type){
          message("trying to convert data into position frequency matrix from count matrix.")
          pcm2pfm(new("pcm", mat=as.matrix(m), name=.f))
        }else{
          stop("Can not import data into PFM. Columns of PFM must add up to 1.0")
        }
      }
    }
  }, m, isCnt, names(m))
}

## output named list of matrix
importFASTAlikeFile <- function(fn, comment.char=">"){
  stopifnot(length(fn)==1)
  lines <- readLines(fn)
  lines <- lines[!grepl("^\\s*$", lines)]
  sep <- grepl(paste0("^", comment.char), lines)
  if(sum(sep)<1){
    ## not fasta like file
    ## insert the filename as fasta like file
    ## check point: only four rows
    if(length(lines)!=4){
      stop("Trying to convert the file into fasta like file.",
           "But the lines in file is not 4.")
    }
    lines <- c(paste0(comment.char, sub("\\.(pcm|pfm|cm|jaspar|scpd|beeml|cisbp)$", "", basename(fn))), lines)
    sep <- grepl(paste0("^", comment.char), lines)
  }
  tfNames <- make.names(sub(paste0("^", comment.char, "\\s*"), "", lines[sep]), unique = TRUE, allow_ = TRUE)
  sep.f <- diff(c(which(sep), length(lines)+1))
  if(any(sep.f!=5)){
    stop("The file contain unexpect lines.",
         "The file should be fasta-like file.",
         "The first line is start with >",
         "and followed with four lines indicates the A, C, G, T counts")
  }
  lines <- lines[!sep]
  seq.f <- rep(tfNames, each=4)
  tfData <- split(lines, seq.f)
  rown <- c("A", "C", "G", "T")
  lapply(tfData, function(.ele){
    ## check ACGT
    .rown <- substr(.ele, 1, 1)
    if(any(grepl("[a-zA-Z]", .rown))){
      if(all(.rown==rown)){
        ## remove the characters
        .ele <- gsub("^[^0-9\\-\\.]*([0-9\\-\\.].*[0-9a-eA-E\\.])[^0-9a-eA-E\\.]*$", "\\1", gsub("\\t", " ", .ele))
      }else{
        if(all(.rown %in% rown)){
          .ele <- .ele[order(.rown)]
          .ele <- gsub("^[^0-9\\-\\.]*([0-9\\-\\.].*[0-9a-eA-E\\.])[^0-9a-eA-E\\.]*$", "\\1", gsub("\\t", " ", .ele))
        }else{
          stop("unexpect rows. The rows should be A, C, G and T.")
        }
      }
    }
    .ele <- matrix(scan(text=.ele, what=double()), nrow=4, byrow = TRUE)
    rownames(.ele) <- rown
    .ele
  })
}
importMulFASTAlike <- function(fns, comment.char=">"){
  d <- lapply(fns, function(.ele){
    importFASTAlikeFile(.ele, comment.char)
  })
  do.call(c, d)
}

importTRANSFAClikeFile <- function(fn){
  stopifnot(length(fn)==1)
  lines <- readLines(fn)
  lines <- lines[!grepl("^\\s*$", lines)]
  ## check PO and XX
  POs <- grepl("^P[O0]\\s", lines)
  XXs <- grepl("^XX", lines)
  if(sum(XXs)<1 || sum(POs)<1){
    ## not TRANSFAC like file
    stop("File is not TRANSFAC like file.")
  }
  ## find the data blocks
  sep <- data.frame(id=c(which(POs), which(XXs)), 
                    group=rep(c("PO","XX"), c(sum(POs), sum(XXs))),
                    stringsAsFactors = FALSE)
  sep <- sep[order(sep$id), ]
  sep <- sep[sep$group=="XX" & c("XX", sep$group[-nrow(sep)])=="PO", "id"]
  sep <- rep(sep, diff(c(0, sep)))
  lines <- split(lines[seq_along(sep)], sep)
  m <- lapply(lines, function(.ele){
    ## check lines with ID, AC, DE, BF
    getMarker <- function(s) sub(paste0("^", s, "\\s+"), "", .ele[grepl(paste0("^", s, "\\s"), .ele)])
    ID <- getMarker("ID")
    AC <- getMarker("AC")
    DE <- getMarker("DE")
    BF <- getMarker("BF")
    tfName <- c(ID, AC, DE, BF)[1]
    ## check P0
    Po <- which(grepl("^P[O0]\\s", .ele))
    .ele <- .ele[Po:(length(.ele)-1)]
    ## check A C G T
    if(!grepl("^P[O0]\\s+A\\s+C\\s+G\\s+T", .ele[1], ignore.case = TRUE)){
      stop("PO line must be in the order of A, C, G and T")
    }
    .ele <- .ele[-1]
    ## remove last column if there is consensus
    .ele <- gsub("[a-zA-Z]*\\s*$", "", .ele)
    .ele <- matrix(scan(text=.ele, what = double()), ncol = 5, byrow=TRUE)
    ## check first column
    if(!all(.ele[, 1]==seq_along(.ele[, 1]))){
      stop("In counts section, rownname should be from 01 to n.")
    }
    .ele <- .ele[, -1]
    colnames(.ele) <- c("A", "C", "G", "T")
    .ele <- t(.ele)
    list(name=tfName, mat=.ele)
  })
  names(m) <- make.names(sapply(m, function(.ele) .ele$name), unique = TRUE, allow_ = TRUE)
  lapply(m, function(.ele) .ele$mat)
}

importMulTRANSFAClike <- function(fns){
  d <- lapply(fns, function(.ele){
    importTRANSFAClikeFile(.ele)
  })
  do.call(c, d)
}

importCisBP <- function(fn){
  stopifnot(length(fn)==1)
  lines <- readLines(fn)
  lines <- lines[!grepl("^\\s*$", lines)]
  ## check Pos and XX
  POs <- grepl("^Pos\\s", lines)
  TFs <- grepl("^TF", lines) & !grepl("Name", lines)
  if(sum(TFs)<1 || sum(POs)<1){
    stop("File is not CisBP PWM file.")
  }
  ## find the data blocks
  sep <- which(TFs)
  if(sep[1]!=1){
    stop("File is not CisBP PWM file.")
  }
  sep <- c(sep, length(lines)+1)
  sep <- diff(sep)
  sep <- rep(seq_along(sep), sep)
  lines <- split(lines[seq_along(sep)], sep)
  m <- lapply(lines, function(.ele){
    ## check lines with ID, AC, DE, BF
    getMarker <- function(s) sub(paste0("^", s, "\\s+"), "", .ele[grepl(paste0("^", s, "\\s"), .ele)])
    TF <- getMarker("TF")[1]
    TFn <- getMarker("TF\\s+Name")
    Motif <- getMarker("Motif")
    tfName <- c(TFn, TF, Motif)[1]
    ## check P0
    Po <- which(grepl("^Pos\\s", .ele))
    .ele <- .ele[Po:(length(.ele)-1)]
    ## check A C G T
    if(!grepl("^Pos\\s+A\\s+C\\s+G\\s+T", .ele[1], ignore.case = TRUE)){
      stop("Pos line must be in the order of A, C, G and T")
    }
    .ele <- .ele[-1]
    ## remove last column if there is consensus
    .ele <- gsub("[a-zA-Z]*\\s*$", "", .ele)
    .ele <- matrix(scan(text=.ele, what = double()), ncol = 5, byrow=TRUE)
    ## check first column
    if(!all(.ele[, 1]==seq_along(.ele[, 1]))){
      stop("In counts section, rownname should be from 01 to n.")
    }
    .ele <- .ele[, -1]
    colnames(.ele) <- c("A", "C", "G", "T")
    .ele <- t(.ele)
    list(name=tfName, mat=.ele)
  })
  names(m) <- make.names(sapply(m, function(.ele) .ele$name), unique = TRUE, allow_ = TRUE)
  lapply(m, function(.ele) .ele$mat)
}

importMulCisBP <- function(fns){
  d <- lapply(fns, function(.ele){
    importCisBP(.ele)
  })
  do.call(c, d)
}


importM_pfm <- function(fns){
  importMulFASTAlike(fns)
}

importM_jaspar <- function(fns){
  importMulFASTAlike(fns)
}

importM_cm <- function(fns){
  importMulFASTAlike(fns)
}

importM_pcm <- function(fns){
  importMulFASTAlike(fns)
}

importM_transfac <- function(fns){
  importMulTRANSFAClike(fns)
}

importM_cisbp <- function(fns){
  importMulCisBP(fns)
}

importM_beeml <- function(fns){
  importMulFASTAlike(fns, comment.char="#")
}

importM_scpd <- function(fns){
  importMulFASTAlike(fns)
}

importM_meme <- function(fns){
  m <- lapply(fns, function(fn){
    lines <- readLines(fn)
    ## MEME version
    if(!grepl("MEME\\s+version\\s", lines[1])){
      stop("MEME version is required in first line.")
    }
    lines <- lines[-1]
    alphabet <- lines[grepl("^ALPHABET=", lines)]
    alphabet <- sub("^.*(AC.*)\\s*$", "\\1", alphabet)
    alphabet <- switch(alphabet, ACGT="DNA", ACGU="RNA", ACDEFGHIKLMNPQRSTVWY="AA", "others")
    letters <- switch(alphabet, DNA=c("A", "C", "G", "T"),
                      RNA=c("A", "C", "G", "U"),
                      AA=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"),
                      others={
                        fr <- which(grepl("^ALPHABET=", lines))
                        to <- which(grepl("^END\\s+ALPHABET", lines))
                        substr(lines[(fr+1):(to-1)], 1, 1)
                      })
    
    tfNames <- sub("^MOTIF\\s+", "", lines[grepl("^MOTIF\\s+", lines)])
    mat.frs <- which(grepl("^letter-probability matrix", lines))
    if(length(tfNames)!=length(mat.frs)){
      stop("Can not convert the MEME file.", fn)
    }
    motifs <- lapply(mat.frs, function(mat.fr){
      mat.info <- lines[mat.fr]
      alength <- sub("^.*?alength=\\s*(\\d+)\\s*w=.*$", "\\1", mat.info)
      w <- sub("^.*?w=\\s*(\\d+)\\s.*$", "\\1", mat.info)
      mat.lines <- lines[(mat.fr+1):(mat.fr+as.numeric(w))]
      mat <- matrix(scan(text=mat.lines, what = double()), ncol = as.numeric(alength), byrow=TRUE)
      mat <- t(mat)
      if(length(letters)!=nrow(mat)){
        stop("alength is not as same as length of letters")
      }
      rownames(mat) <- letters
      mat
    })
    names(motifs) <- make.names(tfNames, unique = TRUE, allow_ = TRUE)
    
    ## background
    background <- NA
    if(any(grepl("^Background\\s+letter\\s+frequencies", lines))){
      bcklines.fr <- which(grepl("^Background\\s+letter\\s+frequencies", lines))+1
      bcklines <- lines[bcklines.fr:(length(lines))]
      bcklines.keep <- which(grepl("^[a-zA-Z]\\s+[0-9\\.]+", bcklines))
      bcklines.keep.diff <- diff(c(0, bcklines.keep))
      bcklines.keep.to <- which(bcklines.keep.diff!=1)
      if(length(bcklines.keep.to)<1){
        bcklines.keep.to <- length(bcklines.keep.to)
      }else{
        bcklines.keep.to <- bcklines.keep.to[1]-1
      }
      bcklines <- bcklines[seq.int(bcklines.keep.to)]
      if(length(bcklines)>0){
        bcklines <- paste(bcklines, collapse = " ")
        bcklines <- sub("\\s*$", "", bcklines)
        bcklines <- strsplit(bcklines, "\\s+")[[1]]
        bcklines <- matrix(bcklines, nrow=2)
        background <- as.numeric(bcklines[2, ])
        names(background) <- bcklines[1, ]
      }
      mapply(function(mat, tfName) new("pfm", mat=mat, name=tfName, alphabet=alphabet, background=background),
             motifs, names(motifs))
    }else{
      mapply(function(mat, tfName) new("pfm", mat=mat, name=tfName, alphabet=alphabet),
             motifs, names(motifs))
    }
  })
  m <- unlist(m)
  m
}

importM_psam <- function(fns){
  m <- lapply(fns, function(fn){
    lines <- xmlToList(xmlParse(fn))
    psam <- trimws(lines$psam) ## now ACGT only, don't know what is exactly format of MatrixREDUCE 2.0
    psam <- read.delim(text = psam, header = FALSE, skip = 2)
    psam <- t(psam)
    rownames(psam) <- c("A", "C", "G", "T")
    new("psam", mat=psam, name=make.names(basename(fn)))
  })
  names(m) <- make.names(sapply(m, function(.ele) .ele@name), unique = TRUE, allow_ = TRUE)
  m
}
