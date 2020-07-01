#' Class \code{"pssm"}
#' 
#' An object of class \code{"pssm"} represents the position specific score 
#' matrix of a DNA/RNA/amino-acid sequence motif. The entry stores a matrix, 
#' which in row i, column j gives the log-odds probability of nucleotide/or 
#' amino acid i in position j of the motif.
#' 
#' 
#' @name pssm-class
#' @aliases pssm
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("pssm", mat, name, alphabet, color, background)}.
#' @keywords classes
#' @export
#' @examples
#' 
#' pcm <- read.table(file.path(find.package("motifStack"), 
#'                             "extdata", "bin_SOLEXA.pcm"))
#' pcm <- pcm[,3:ncol(pcm)]
#' rownames(pcm) <- c("A","C","G","T")
#' motif <- pcm2pssm(pcm)
#' motif <- new("pssm", mat=motif, name="bin_SOLEXA")
#' plot(motif)
#' 
setClass("pssm", 
         representation(mat="matrix", name="character", alphabet="character", 
                        color="character", background="numeric", tags="list",
                        markers="list"),
         validity=function(object){
           re<-TRUE
           if (is.null(rownames(object@mat))) re<-"rownames of pssm is empty"
           if (is.null(names(object@color))) re<-"names of color is empty"
           if (length(setdif(rownames(object@mat), names(object@color)))>0) re<-"not every alphabet in pssm has a corresponding color"
           if (is.null(names(object@background))) re<-"names of background is empty"
           if (length(setdif(rownames(object@mat), names(object@background)))>0) re<-"not every alphabet in pssm has a background"
           if (!object@alphabet %in% c("DNA", "RNA", "AA", "others")) re<-"alphabet must be either DNA, RNA, AA or others"
           if (object@alphabet=="DNA" & !all(rownames(object@mat) %in% c("A","C","G","T"))) re<-"rownames of pssm for DNA motifs must be capital A, C, G and T"
           if (object@alphabet=="RNA" & !all(rownames(object@mat) %in% c("A","C","G","U"))) re<-"rownames of pssm for RNA motifs must be capital A, C, G and U"
           if (object@alphabet=="AA" & !all(rownames(object@mat) %in% LETTERS[1:26])) re<-"rownames of pssm for AA motifs must be capital letters"
           if (any(!is(object@markers, "marker"))) re<- "markers must be a list of marker object"
           re
         }
)

#' "pssm" methods
#' 
#' methods for pssm objects.
#' 
#' 
#' @rdname pssm-class
#' @aliases $,pssm-method $<-,pssm-method
#' @docType methods
#' @param x An object of class \code{pssm}. For \code{getIC}, if parameter p is
#' followed, x should be an object of matrix.
#' @param name Slot name.
#' @param y Not use.
#' @param n how many spaces should be added.
#' @param b logical value to indicate where the space should be added.
#' @param \dots Further potential arguments passed to \code{plotMotifLogo}.
#' @param row.names,optional see as.data.frame
#' @section Methods: \describe{ \item{addBlank}{\code{signature(x="pssm",
#' n="numeric", b="logical")} add space into the position frequency matrix for
#' alignment. b is a bool value, if TRUE, add space to the 3' end, else add
#' space to the 5' end. n indicates how many spaces should be added.}
#' 
#' \item{matrixReverseComplement}{\code{signature(x = "pssm")} get the reverse
#' complement of position frequency matrix.}
#' 
#' \item{plot}{\code{signature(x = "pssm")} Plots the sequence logo of the
#' position frequency matrix. }
#' 
#' \item{$, $<-}{Get or set the slot of \code{\link{pssm-class}}}
#' \item{as.data.frame}{convert \code{\link{pssm-class}} to a data.frame}
#' \item{format}{return the name_pssm of \code{\link{pssm-class}}} }
#' @keywords classes
#' @examples
#' 
#' pcm <- read.table(file.path(find.package("motifStack"), 
#'                   "extdata", "bin_SOLEXA.pcm"))
#' pcm <- pcm[,3:ncol(pcm)]
#' rownames(pcm) <- c("A","C","G","T")
#' motif <- pcm2pssm(pcm)
#' motif <- new("pssm", mat=motif, name="bin_SOLEXA")
#' matrixReverseComplement(motif)
#' addBlank(motif, 1, FALSE)
#' addBlank(motif, 3, TRUE)
#' as(motif,"matrix")
#' as.data.frame(motif)
#' format(motif)
#' 
#' @exportMethod `$` `$<-`
setMethod("$", "pssm", function(x, name) slot(x, name))
setReplaceMethod("$", "pssm",
                 function(x, name, value){
                   slot(x, name, check = TRUE) <- value
                   x
                 })
#' @importFrom grDevices rainbow
setMethod("initialize","pssm",function(.Object, mat, name, alphabet, 
                                       color, background, tags, markers){
  if(mode(mat)!="numeric") stop("mat must be a numeric matrix")
  if(is.null(rownames(mat))) stop("rownames of pssm is empty")
  rownames(mat)<-toupper(rownames(mat))
  rname<-rownames(mat)
  if(missing(alphabet)){
    rname<-rname[order(rname)]
    if(all(rname %in% c("A","C","G","T"))) alphabet<-"DNA"
    else{
      if(all(rname %in% c("A","C","G","U"))) alphabet<-"RNA"
      else{
        if(all(rname %in% LETTERS[1:26]) & nrow(mat)==20) alphabet<-"AA"
        else alphabet<-"others"
      }
    }
  }
  if(missing(color)){
    if(alphabet!="others"){
      color<-colorset(alphabet, "auto")
    }else{
      color <- rainbow(nrow(mat))
      names(color) <- rname
    }
  }
  names(color) <- toupper(names(color))
  if(missing(background)){
    background<-rep(1/nrow(mat), nrow(mat))
    names(background)<-rname
  }
  if(length(background)<nrow(mat)){
    background <- rep(background, nrow(mat))[1:nrow(mat)]
    names(background) <- rname
  }
  .Object@mat            =    mat
  .Object@name        =    name
  .Object@alphabet    =    alphabet
  .Object@color        =    color
  .Object@background    =    background
  if(!missing(tags)) .Object@tags = tags
  if(!missing(markers)) .Object@markers = markers
  .Object
})

#' @rdname pssm-class 
#' @aliases plot,pssm,ANY-method
#' @exportMethod plot
setMethod("plot", signature(x="pssm"), 
          function(x, y="missing", ...){
            plotAffinityLogo(psam=x, ...)
          }
)

#' @name coerce()
#' @rdname pssm-class
#' @aliases coerce,pssm,matrix-method
#' @exportMethod coerce
setAs(from="pssm", to="matrix", function(from){
  from@mat
})


#' @rdname pssm-class
#' @exportMethod matrixReverseComplement
#' @aliases matrixReverseComplement,pssm-method
setMethod("matrixReverseComplement", "pssm", function(x){
  if(x@alphabet!="DNA") stop("alphabet of pssm must be DNA")
  mat<-x@mat
  ## double check the rownames to A, C, G, T.
  if(!all(toupper(rownames(mat))==c("A", "C", "G", "T"))){
    mat <- mat[match(c("A", "C", "G", "T"), toupper(rownames(mat))), ]
  }
  rc<-matrix(nrow=nrow(mat),ncol=ncol(mat))
  complements<-c(4,3,2,1)
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      rc[complements[i],(ncol(mat)+1-j)]<-mat[i,j]
    }
  }
  rownames(rc)<-rownames(mat)
  x@mat<-rc
  x
})

#' @rdname pssm-class
#' @exportMethod addBlank
#' @aliases addBlank,pssm,numeric,logical-method
setMethod("addBlank", signature(x="pssm", n="numeric", b="logical"), function(x, n, b){
  if(x@alphabet!="DNA") stop("alphabet of pssm must be DNA")
  N<-matrix(rep(0, n*4), nrow=4)
  if(b){
    N<-cbind(x@mat,N)
  }else{
    N<-cbind(N,x@mat)
    if(length(x@markers)>0){
      x@markers <- lapply(x@markers, function(.ele){
        .ele@start <- .ele@start+n
        .ele@stop <- .ele@stop+n
        .ele
      })
    }
  }
  rownames(N)<-rownames(x@mat)
  x@mat<-N
  x
})

## for data.frame
#' @rdname pssm-class
#' @exportMethod as.data.frame
#' @aliases as.data.frame,pssm-method
setMethod("as.data.frame", signature(x="pssm"), function(x, row.names = NULL, optional = FALSE, ...){
  as.data.frame(x@mat, ...)
})

#' @rdname pssm-class
#' @exportMethod format
#' @aliases format,pssm-method
setMethod("format", signature(x="pssm"), function(x, ...){
  paste0(x@name, "_pssm")
})
