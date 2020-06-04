#' Class \code{marker}
#' 
#' An object of class \code{"marker"} represents a marker in a motif
#' 
#' 
#' @name marker-class
#' @aliases marker marker-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("marker", type, start, stop, label, gp)}.
#' @keywords classes
#' @exportClass marker
#' @examples
#' 
#'  new("marker", type="rect", start=c(2, 4), gp=gpar(lty=3))
#' 
# Import some S3 classes
setOldClass("gpar")

setClass(
  "marker", 
  representation(type="character", start="numeric", 
                 stop="numeric", label="character", 
                 gp="gpar"),
  validity=function(object){
    re<-TRUE
    if (length(object@type)!=1) 
      re<-"type should be 1 of 'rect', 'text', 'line'"
    if (!object@type%in%c('rect', 'text', 'symbol', 'line')) 
      re<-"type should be 1 of 'rect', 'text', 'line'"
    if (any(object@stop<object@start)) 
      re <- "stop must be no less than start"
    if (any(start!=round(start, digits = 0)) || 
        any(stop!=round(stop, digits = 0))){
      re <- "start and stop must be an integer."
    }
    re
  }
)

#' @rdname marker-class
#' @docType methods
#' @param x A marker object
#' @param name slot name of marker object
#' @exportMethod `$` `$<-`
#' @aliases $,marker-method $<-,marker-method
setMethod("$", "marker", function(x, name) slot(x, name))
setReplaceMethod("$", "marker",
                 function(x, name, value){
                   slot(x, name, check = TRUE) <- value
                   x
                 })

#' @importFrom grid gpar
setMethod("initialize","marker",function(.Object, type, start, stop, label, gp){
  if (length(type)!=1) 
    re<-"type should be 1 of 'rect', 'text', 'line'"
  if (!type%in%c('rect', 'text', 'symbol', 'line'))
    re<-"type should be 1 of 'rect', 'text', 'line'"
  if(missing(label)) label = character()
  if(missing(start)) start = numeric()
  if(missing(stop)) stop = start
  if(missing(gp)) gp = gpar()
  .Object@type            =    type;
  .Object@start        =    start;
  .Object@stop      = stop;
  .Object@label    =    label;
  .Object@gp        =    gp;
  .Object
})

