#' Class \code{ouNode}
#' 
#' An object of class \code{"ouNode"} represents a motif node in a cluster tree
#' 
#' 
#' @name ouNode-class
#' @aliases ouNode ouNode-class $,ouNode-method $<-,ouNode-method
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("ouNode", left, right, parent, distl, distr, sizel, sizer)}.
#' @keywords classes
#' @export
#' @import methods
#' @examples
#' 
#'  new("ouNode", left="A", right="B", parent="Root", distl=1, distr=2, sizel=1, sizer=1)
#' 

setClass("ouNode", representation(left="character", 
                                  right="character", 
                                  parent="character", 
                                  distl="numeric", 
                                  distr="numeric",
                                  sizel="numeric",
                                  sizer="numeric"))
#' @rdname ouNode-class
#' @docType methods
#' @param x A ouNode object
#' @param name slot name of ouNode object
#' @exportMethod `$` `$<-`
setMethod("$", "ouNode", function(x, name) slot(x, name))
setReplaceMethod("$", "ouNode",
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })

#' Class \code{"motifSig"}
#' 
#' An object of class \code{"motifSig"} represents the output of function
#' \link{motifSignature}
#' 
#' 
#' @name motifSig-class
#' @aliases motifSig
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("motifSig", signature, freq, nodelist, gpcol)}.
#' @keywords classes
#' @export
setClass("motifSig", representation(signatures="list", freq="numeric", 
                                    nodelist="list", gpcol="character"),
         validity=function(object){
             re<-TRUE
             if (any(unlist(lapply(object@nodelist, function(.ele) !inherits(.ele, "ouNode"))))) re<-"nodelist must be a list of ouNode objects"
             if (any(unlist(lapply(object@signatures, function(.ele) !inherits(.ele, "pfm"))))) re<-"signatures must be a list of pfm objects"
             if (length(object@signatures) != length(object@freq)) re<-"ength of signatures should be same as freq"
             re
         }
)

#' "motifSig" methods
#' 
#' methods for motifSig objects.
#' 
#' 
#' @rdname motifSig-class
#' @aliases signatures signatures,motifSig-method frequence
#' @docType methods
#' @param object An object of class \code{motifSig}.
#' @section Methods: \describe{ \item{signatures}{\code{signature(object =
#' "motifSig")} return the signatures of motifSig }
#' 
#' \item{frequence}{\code{signature(object = "motifSig")} return the frequency
#' of motifSig}
#' 
#' \item{nodelist}{\code{signature(object = "motifSig")} return the nodelist of
#' motifSig}
#' 
#' \item{sigColor}{\code{signature(object = "motifSig")} return the group color
#' sets of motifSig}
#' 
#' \item{$, $<-}{Get or set the slot of \code{\link{motifSig}}} }
#' @keywords classes
#' @exportMethod signatures
setGeneric("signatures", function(object) standarGeneric("signatures"))
setMethod("signatures", signature(object="motifSig"), function(object) object@signatures)

#' @rdname motifSig-class
#' @exportMethod frequence
#' @aliases frequence,motifSig-method
setGeneric("frequence", function(object) standarGeneric("frequence"))
setMethod("frequence", signature(object="motifSig"), function(object) object@freq)

#' @rdname motifSig-class
#' @exportMethod nodelist
#' @aliases nodelist nodelist,motifSig-method
setGeneric("nodelist", function(object) standarGeneric("nodelist"))
setMethod("nodelist", signature(object="motifSig"), function(object) object@nodelist)

#' @rdname motifSig-class
#' @exportMethod sigColor
#' @aliases sigColor sigColor,motifSig-method
setGeneric("sigColor", function(object) standarGeneric("sigColor"))
setMethod("sigColor", signature(object="motifSig"), function(object) object@gpcol)

#' @rdname motifSig-class
#' @exportMethod `$` `$<-`
#' @param x A motifSig object
#' @param name slot name of motifSig object
#' @aliases $,motifSig-method $<-,motifSig-method
setMethod("$", "motifSig", function(x, name) slot(x, name))
setReplaceMethod("$", "motifSig",
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })
