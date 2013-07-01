checkInteger <- function(N){
    !length(grep("[^[:digit:]]", as.character(N)))
}

hex2psrgb<-function(col){
    col<-grDevices::col2rgb(col)
    col<-col/255
    col<-paste(col,collapse=" ")
    col
}

highlightCol <- function (col, alpha=0.5){
	n <- names(col)
	col <- grDevices::col2rgb(col,alpha=T)
	col <- apply(col, 2, function(.ele, alpha){rgb(.ele[1], .ele[2], .ele[3], alpha=ceiling(alpha*.ele[4]), maxColorValue=255)}, alpha)
	col <- unlist(col)
	names(col) <- n
	col
}

motifStack_private_fontsize <- 72
coloredSymbols <- function(ncha, font, color, rname, fontsize=motifStack_private_fontsize){
	symbols<-list()
	for(i in 1:ncha){
		ps<-paste("%!PS\n/",font," findfont\n",fontsize," scalefont\n",
				  motifStack:::hex2psrgb(color[i])," setrgbcolor\nsetfont\nnewpath\n0 0 moveto\n(",
				  rname[i],") show",sep="")
		psfilename<-tempfile()
		psfilename <- gsub("\\", "/", psfilename, fixed=TRUE)
# step1 create a ps file
		cat(ps,file=paste(psfilename,".ps",sep=""))
# step2 convert it by grImport::PostScriptTrace
		grImport::PostScriptTrace(paste(psfilename,".ps",sep=""), paste(psfilename,".xml",sep=""))
# step3 read by grImport::readPicture
		symbols[[i]]<-grImport::readPicture(paste(psfilename,".xml",sep=""))
		unlink(c(paste(psfilename,".ps",sep=""), 
				 paste("capture",basename(psfilename),".ps",sep=""), 
				 paste(psfilename,".xml",sep="")))
	}
	symbols
}

addPseudolog2<-function(x){
    ifelse(x==0, -10, log2(x))
}

## get Information Entropy from matrix
getIE<-function(x){
    addPseudolog2(nrow(x))
}

UngappedAlignment<-function(pfms, i, threshold, minimalConsensus=0, rcpostfix="(RC)", revcomp=TRUE){
	if(class(pfms[[i]])!="pfm"){
	   pcms <- pfms
	   pfms <- lapply(pfms, pcm2pfm)
	}else{
	   pcms <- NULL
	}
    res<-getAlignedICWithoutGap(pfms[[i-1]], pfms[[i]], threshold, revcomp)
	if(!is.null(pcms)) pfms <- pcms
	if(res$max>=minimalConsensus){
		if(res$rev){
			pfms[[i]]<-matrixReverseComplement(pfms[[i]])
			pfms[[i]]@name<-paste(pfms[[i]]@name, rcpostfix, sep="")
		}
		if(res$offset>0){
			pfms[[i]]<-addBlank(pfms[[i]], res$offset, FALSE)
		}else{
			if(res[1]<0){
				pfms[1:(i-1)]<-lapply(pfms[1:(i-1)],function(.ele) addBlank(.ele, -1*res$offset, FALSE))
			}
		}
	}
    pfms
}
getAlignedICWithoutGap<-function(pfm1, pfm2, threshold, revcomp=TRUE){
	if(class(pfm1)!="pfm" | class(pfm2)!="pfm") stop("class of pfm1 and pfm2 must be pfm")
    offset1<-getoffsetPosByIC(pfm1, pfm2, threshold)
    if(revcomp){
        pfm3<-matrixReverseComplement(pfm2)
        offset2<-getoffsetPosByIC(pfm1, pfm3, threshold)
    }else{
        offset2<-list(k=0,max=0)
    }
    offset<-0
    rev<-FALSE
    if(offset1$max < offset2$max){
        rev<-TRUE
		max<-offset2$max
        offset<-offset2$k
    }else{
		max<-offset1$max
        offset<-offset1$k
    }
    list(offset=offset,rev=rev,max=max)    
}
getICbyBase<-function(p, pfmj){
    re<-rep(0, 5)
    re[1:4]<-pfmj * (addPseudolog2(pfmj) - addPseudolog2(p))
    re[5]<-sum(re[1:4])
    re
}
getALLRbyBase <- function(b, pfmi, pfmj){##return 2x ALLR
    sum(pfmj * (addPseudolog2(pfmi) - addPseudolog2(b))) 
            + sum(pfmi * (addPseudolog2(pfmj) - addPseudolog2(b)))
}
getoffsetPosByIC<-function(pfm1, pfm2, threshold){
	if(class(pfm1)!="pfm" | class(pfm2)!="pfm") stop("class of pfm1 and pfm2 must be pfm")
    score1<-rep(0, ncol(pfm1@mat))
    score2<-rep(0, ncol(pfm2@mat))
    for(i in 1:ncol(pfm1@mat)){
        J<-ifelse(ncol(pfm1@mat)+1-i>ncol(pfm2@mat),ncol(pfm2@mat),ncol(pfm1@mat)+1-i)
        for(j in 1:J){
            ic1<-getICbyBase(pfm1@background, pfm1@mat[,i+j-1])
            ic2<-getICbyBase(pfm2@background, pfm2@mat[,j])
            ic3<-getALLRbyBase(pfm1@background, pfm1@mat[,i+j-1], pfm2@mat[,j])
            if(ic1[5]>threshold & ic2[5]>threshold){
 #               a<-ic1[1:4]>mean(ic1[1:4])
 #               b<-ic2[1:4]>mean(ic2[1:4])
 #               if(any((a&b))) score1[i]<-score1[i]+1
                if(ic3>0) score1[i]<-score1[i]+1
            }
        }
    }
    for(i in 1:ncol(pfm2@mat)){
        J<-ifelse(ncol(pfm2@mat)+1-i>ncol(pfm1@mat),ncol(pfm1@mat),ncol(pfm2@mat)+1-i)
        for(j in 1:J){
            ic2<-getICbyBase(pfm2@background, pfm2@mat[,i+j-1])
            ic1<-getICbyBase(pfm1@background, pfm1@mat[,j])
            ic3<-getALLRbyBase(pfm1@background, pfm1@mat[,j], pfm2@mat[,i+j-1])
            if(ic1[5]>threshold & ic2[5]>threshold){
 #               a<-ic1[1:4]>mean(ic1[1:4])
 #               b<-ic2[1:4]>mean(ic2[1:4])
 #               if(any((a&b))) score2[i]<-score2[i]+1
                if(ic3>0) score2[i]<-score2[i]+1
            }
        }
    }
    k1<-match(max(score1),score1)
    k2<-match(max(score2),score2)
    k=ifelse(score1[k1]<score2[k2], -1*(k2-1), (k1-1))
    max=ifelse(score1[k1]<score2[k2], score2[k2], score1[k1])
    list(k=k,max=max)
}

setClass("Rect", 
		representation(x="numeric", y="numeric", width="numeric", height="numeric"), 
		prototype=prototype(x=0, y=0, width=0, height=0)
)

setGeneric("isContainedIn", function(a, b) standardGeneric("isContainedIn"))

setMethod("isContainedIn", signature(a="Rect", b="Rect"), function(a, b){
	a@x >= b@x && a@y >= b@y && a@x+a@width <= b@x+b@width && a@y+a@height <= b@y+b@height
})

setClass("pos",
	representation(box="Rect", beta="numeric", sig="pfm", freq="numeric", norm="numeric")
)

getPFMid <- function(pfms, nodename, rcpostfix="(RC)"){
	pfmNames <- as.character(unlist(lapply(pfms, function(.ele) .ele@name)))
	pfmNames <- gsub(rcpostfix,"",pfmNames, fixed=T)
	which(pfmNames==nodename)
}

getParentNode <- function(nodelist, nodename){
	for(i in 1:length(nodelist)){
		currNode <- nodelist[[i]]
		if(currNode@left==nodename) return(c(currNode@parent, "left"))
		if(currNode@right==nodename) return(c(currNode@parent, "right"))
	}
	NULL
}