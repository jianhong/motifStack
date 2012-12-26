setClass("pcm", 
		representation(mat="matrix", name="character", alphabet="character", color="character", background="numeric"),
		validity=function(object){
			re<-TRUE
			if (checkInteger(unlist(object@mat))) re<-"Columns of pcm must add up to 1.0"
			if (is.null(rownames(object@mat))) re<-"rownames of pcm is empty"
			if (is.null(names(object@color))) re<-"names of color is empty"
			if (length(setdif(rownames(object@mat), names(object@color)))>0) re<-"not every alphabet in pcm has a corresponding color"
			if (is.null(names(object@background))) re<-"names of background is empty"
			if (length(setdif(rownames(object@mat), names(object@background)))>0) re<-"not every alphabet in pcm has a background"
			if (!object@alphabet %in% c("DNA", "RNA", "AA", "others")) re<-"alphabet must be either DNA, RNA, AA or others"
			if (object@alphabet=="DNA" & !all(rownames(object@mat) %in% c("A","C","G","T"))) re<-"rownames of pcm for DNA motifs must be capital A, C, G and T"
			if (object@alphabet=="RNA" & !all(rownames(object@mat) %in% c("A","C","G","U"))) re<-"rownames of pcm for RNA motifs must be capital A, C, G and U"
			if (object@alphabet=="AA" & !all(rownames(object@mat) %in% LETTERS[1:26])) re<-"rownames of pcm for AA motifs must be capital letters"
			re
		}
)

setMethod("initialize","pcm",function(.Object, mat, name, alphabet, color, background){
	if(is.null(rownames(mat))) stop("rownames of pcm is empty")
	rownames(mat)<-toupper(rownames(mat))
	if(missing(alphabet)){
		rname<-rownames(mat)
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
	if(missing(background)){
		background<-rep(1/nrow(mat), nrow(mat))
		names(background)<-rname
	}
	.Object@mat			=	mat;
	.Object@name		=	name;
	.Object@alphabet	=	alphabet;
	.Object@color		=	color;
	.Object@background	=	background;
	.Object
})

setMethod("plot", signature(x="pcm"), 
	function(x, y="missing", ...){
		plotMotifLogo(pfm=pcm2pfm(x@mat), 
					  motifName=x@name, 
					  p=x@background[rownames(x@mat)], 
					  colset=x@color[rownames(x@mat)], 
					  ...)
	}
)

setAs(from="pcm", to="matrix", function(from){
	from@mat
})


setGeneric("addBlank", function(x, n, b) standardGeneric("addBlank"))
setGeneric("getIC", function(x, p) standardGeneric("getIC"))
setGeneric("matrixReverseComplement", function(x) standardGeneric("matrixReverseComplement"))

setMethod("matrixReverseComplement", "pcm", function(x){
	if(x@alphabet!="DNA") stop("alphabet of pcm must be DNA")
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

setMethod("addBlank", signature(x="pcm", n="numeric", b="logical"), function(x, n, b){
	if(x@alphabet!="DNA") stop("alphabet of pfm must be DNA")
	N<-matrix(rep(rep(0,4), n), nrow=4)
	if(b){
        N<-cbind(x@mat,N)
    }else{
        N<-cbind(N,x@mat)
    }
    rownames(N)<-rownames(x@mat)
    x@mat<-N
	x
})

## get information content profile from pcm
setMethod("getIC", signature(x="pcm"), function(x, p="missing"){
	mat<-pcm2pfm(x@mat)
	colSums(mat * (addPseudolog2(mat) - addPseudolog2(x@background)))
})

setGeneric("pcm2pfm", function(x, background) standardGeneric("pcm2pfm"))

## pcm to pfm from matrix
setMethod("pcm2pfm", signature(x="matrix", background="numeric"), function(x, background){
	s<-apply(x,2,sum)
	s[s==0] <- NA
	p<-apply(x,1,function(.ele) .ele/s)
	n<-length(s[is.na(s)])
	if(any(is.na(s))) p[is.na(p[,1]),] <- c(rep(background, n))
	t(p)
})

setMethod("pcm2pfm", signature(x="matrix"), function(x, background="missing"){
	pcm2pfm(x, rep(1/nrow(x), nrow(x)))
})

setMethod("pcm2pfm", signature(x="data.frame", background="numeric"), function(x, background){
	pcm2pfm(as.matrix(x), background)
})

setMethod("pcm2pfm", signature(x="data.frame"), function(x, background="missing"){
	pcm2pfm(as.matrix(x), rep(1/nrow(x), nrow(x)))
})

## pcm to pfm from pcm object
setMethod("pcm2pfm", signature(x="pcm"), function(x, background="missing"){
	.mat <- pcm2pfm(x@mat, x@background)
	new("pfm", mat=.mat, name=x@name, alphabet=x@alphabet, color=x@color, background=x@background)
})