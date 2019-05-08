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

setMethod("plot", signature(x="psam", y="ANY"), 
          function(x, y="missing", ...){
            plotAffinityLogo(psam=x, ...)
          }
)

setAs(from="psam", to="matrix", function(from){
  from@mat
})


setMethod("matrixReverseComplement", "psam", function(x){
  if(x@alphabet!="DNA") stop("alphabet of psam must be DNA")
  mat<-x@mat
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

setMethod("addBlank", signature(x="psam", n="numeric", b="logical"), function(x, n, b){
  if(x@alphabet!="DNA") stop("alphabet of psam must be DNA")
  N<-matrix(rep(0, n*4), nrow=4)
  if(b){
    N<-cbind(x@mat,N)
  }else{
    N<-cbind(N,x@mat)
  }
  rownames(N)<-rownames(x@mat)
  x@mat<-N
  x
})

## for data.frame
setMethod("as.data.frame", signature(x="psam"), function(x, row.names = NULL, optional = FALSE, ...){
  as.data.frame(x@mat, ...)
})

setMethod("format", signature(x="psam"), function(x, ...){
  paste0(x@name, "_psam")
})
