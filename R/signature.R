setClass("ouNode", representation(left="character", 
                                  right="character", 
                                  parent="character", 
                                  distl="numeric", 
                                  distr="numeric",
                                  sizel="numeric",
                                  sizer="numeric"))

setMethod("$", "ouNode", function(x, name) slot(x, name))
setReplaceMethod("$", "ouNode",
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })

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

setGeneric("signatures", function(object) standarGeneric("signatures"))
setMethod("signatures", signature(object="motifSig"), function(object) object@signatures)

setGeneric("frequence", function(object) standarGeneric("frequence"))
setMethod("frequence", signature(object="motifSig"), function(object) object@freq)

setGeneric("nodelist", function(object) standarGeneric("nodelist"))
setMethod("nodelist", signature(object="motifSig"), function(object) object@nodelist)

setGeneric("sigColor", function(object) standarGeneric("sigColor"))
setMethod("sigColor", signature(object="motifSig"), function(object) object@gpcol)

setMethod("$", "motifSig", function(x, name) slot(x, name))
setReplaceMethod("$", "motifSig",
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })
