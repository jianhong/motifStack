#' Class \code{"psam"}
#' 
#' An object of class \code{"psam"} represents the position specific affinity
#' matrix (PSAM) of a DNA/RNA/amino-acid sequence motif. The entry stores a
#' matrix, which in row i, column j gives the affinity of observing
#' nucleotide/or amino acid i in position j of the motif.
#' 
#' 
#' @name psam-class
#' @aliases psam
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("psam", mat, name, alphabet, color)}.
#' @keywords classes
#' @export
#' @examples
#' 
#' motif <- importMatrix(file.path(find.package("motifStack"), "extdata", "PSAM.mxr"), 
#'          format="psam")[[1]]
#' plot(motif)
#' 
setClass("psam", 
         representation(mat="matrix", name="character", alphabet="character", 
                        color="character", tags="list",
                        markers="list"),
         validity=function(object){
           re<-TRUE
           if (abs(1-max(object@mat)) > 0.01) re<-"PSAM matrix should be normalized with respect to the value of the highest affinity oligonucleotied."
           if (any(object@mat<0)) re<-"PSAM matrix contain value less than 0"
           if (is.null(rownames(object@mat))) re<-"rownames of PSAM is empty"
           if (is.null(names(object@color))) re<-"names of color is empty"
           if (length(setdif(rownames(object@mat), names(object@color)))>0) re<-"not every alphabet in psam has a corresponding color"
           if (!object@alphabet %in% c("DNA", "RNA", "AA", "others")) re<-"alphabet must be either DNA, RNA, AA or others"
           if (object@alphabet=="DNA" & !all(rownames(object@mat) %in% c("A","C","G","T"))) re<-"rownames of PSAM for DNA motifs must be capital A, C, G and T"
           if (object@alphabet=="RNA" & !all(rownames(object@mat) %in% c("A","C","G","U"))) re<-"rownames of PSAM for RNA motifs must be capital A, C, G and U"
           if (object@alphabet=="AA" & !all(rownames(object@mat) %in% LETTERS[1:26])) re<-"rownames of PSAM for AA motifs must be capital letters"
           if (any(!is(object@markers, "marker"))) re<- "markers must be a list of marker object"
           re
         }
)

#' "psam" methods
#' 
#' methods for psam objects.
#' 
#' 
#' @name psam-class
#' @aliases $,psam-method $<-,psam-method
#' @docType methods
#' @param x An object of class \code{psam}.
#' @param name Slot name.
#' @param y Not use.
#' @param n how many spaces should be added.
#' @param b logical value to indicate where the space should be added.
#' @param \dots Further potential arguments passed to \code{plotAffinityLogo}.
#' @param row.names,optional see as.data.frame
#' @section Methods: \describe{ \item{addBlank}{\code{signature(x="psam",
#' n="numeric", b="logical")} add space into the position specific affinity
#' matrix for alignment. b is a bool value, if TRUE, add space to the 3' end,
#' else add space to the 5' end. n indicates how many spaces should be added.}
#' 
#' \item{matrixReverseComplement}{\code{signature(x = "psam")} get the reverse
#' complement of position specific affinity matrix.}
#' 
#' \item{plot}{\code{signature(x = "psam")} Plots the affinity logo of the
#' position specific affinity matrix. }
#' 
#' \item{$, $<-}{Get or set the slot of \code{\link{psam-class}}}
#' \item{as.data.frame}{convert \code{\link{psam-class}} to a data.frame}
#' \item{format}{return the name_pfm of \code{\link{psam-class}}} }
#' @keywords classes
#' @examples
#' 
#' motif <- importMatrix(file.path(find.package("motifStack"), "extdata", "PSAM.mxr"), 
#'                       format="psam")[[1]]
#' matrixReverseComplement(motif)
#' addBlank(motif, 1, FALSE)
#' addBlank(motif, 3, TRUE)
#' as(motif,"matrix")
#' as.data.frame(motif)
#' format(motif)
#' 
#' @exportMethod `$` `$<-`
setMethod("$", "psam", function(x, name) slot(x, name))
setReplaceMethod("$", "psam",
                 function(x, name, value){
                   slot(x, name, check = TRUE) <- value
                   x
                 })

setMethod("initialize","psam",function(.Object, mat, name, alphabet, color, tags, markers){
  if(mode(mat)!="numeric") stop("mat must be a numeric matrix")
  if(is.null(rownames(mat))) stop("rownames of PSAM is empty")
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
    color<-colorset(alphabet, "auto")
  }
  .Object@mat            =    mat
  .Object@name        =    name
  .Object@alphabet    =    alphabet
  .Object@color        =    color
  if(!missing(tags)) .Object@tags = tags
  if(!missing(markers)) .Object@markers = markers
  .Object
})

#' @rdname psam-class
#' @aliases plot,psam,ANY-method
#' @exportMethod plot
setMethod("plot", signature(x="psam"), 
          function(x, y="missing", ...){
            plotAffinityLogo(psam=x, ...)
          }
)

#' @name as()
#' @rdname psam-class
#' @aliases coerce,psam,matrix-method
#' @exportMethod coerce
setAs(from="psam", to="matrix", function(from){
  from@mat
})

#' @rdname psam-class
#' @exportMethod matrixReverseComplement
#' @aliases matrixReverseComplement,psam-method
setMethod("matrixReverseComplement", "psam", function(x){
  if(x@alphabet!="DNA") stop("alphabet of psam must be DNA")
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

#' @rdname psam-class
#' @exportMethod addBlank
#' @aliases addBlank,psam,numeric,logical-method
setMethod("addBlank", signature(x="psam", n="numeric", b="logical"), function(x, n, b){
  if(x@alphabet!="DNA") stop("alphabet of psam must be DNA")
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
#' @rdname psam-class
#' @exportMethod as.data.frame
#' @aliases as.data.frame,psam-method
setMethod("as.data.frame", signature(x="psam"), function(x, row.names = NULL, optional = FALSE, ...){
  as.data.frame(x@mat, ...)
})

#' @rdname psam-class
#' @exportMethod format
#' @aliases format,psam-method
setMethod("format", signature(x="psam"), function(x, ...){
  paste0(x@name, "_psam")
})
