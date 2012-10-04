hex2psrgb<-function(col){
    col<-grDevices::col2rgb(col)
    col<-col/255
    col<-paste(col,collapse=" ")
    col
}

addPseudolog2<-function(x){
    ifelse(x==0, -10, log2(x))
}

## get Information Entropy from matrix
getIE<-function(x){
    addPseudolog2(nrow(x))
}

UngappedAlignment<-function(pfms, i, threshold, minimalConsensus=0){
    res<-getAlignedICWithoutGap(pfms[[i-1]], pfms[[i]], threshold)
	if(res$max>=minimalConsensus){
		if(res$rev){
			pfms[[i]]<-matrixReverseComplement(pfms[[i]])
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
getAlignedICWithoutGap<-function(pfm1, pfm2, threshold){
	if(class(pfm1)!="pfm" | class(pfm2)!="pfm") stop("class of pfm1 and pfm2 must be pfm")
    pfm3<-matrixReverseComplement(pfm2)
    offset1<-getoffsetPosByIC(pfm1, pfm2, threshold)
    offset2<-getoffsetPosByIC(pfm1, pfm3, threshold)
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
getoffsetPosByIC<-function(pfm1, pfm2, threshold){
	if(class(pfm1)!="pfm" | class(pfm2)!="pfm") stop("class of pfm1 and pfm2 must be pfm")
    score1<-rep(0, ncol(pfm1@mat))
    score2<-rep(0, ncol(pfm2@mat))
    for(i in 1:ncol(pfm1@mat)){
        J<-ifelse(ncol(pfm1@mat)+1-i>ncol(pfm2@mat),ncol(pfm2@mat),ncol(pfm1@mat)+1-i)
        for(j in 1:J){
            ic1<-getICbyBase(pfm1@background, pfm1@mat[,i+j-1])
            ic2<-getICbyBase(pfm2@background, pfm2@mat[,j])
            if(ic1[5]>threshold & ic2[5]>threshold){
                a<-ic1[1:4]>mean(ic1[1:4])
                b<-ic2[1:4]>mean(ic2[1:4])
                if(any((a&b))) score1[i]<-score1[i]+1
            }
        }
    }
    for(i in 1:ncol(pfm2@mat)){
        J<-ifelse(ncol(pfm2@mat)+1-i>ncol(pfm1@mat),ncol(pfm1@mat),ncol(pfm2@mat)+1-i)
        for(j in 1:J){
            ic2<-getICbyBase(pfm2@background, pfm2@mat[,i+j-1])
            ic1<-getICbyBase(pfm1@background, pfm1@mat[,j])
            if(ic1[5]>threshold & ic2[5]>threshold){
                a<-ic1[1:4]>mean(ic1[1:4])
                b<-ic2[1:4]>mean(ic2[1:4])
                if(any((a&b))) score2[i]<-score2[i]+1
            }
        }
    }
    k1<-match(max(score1),score1)
    k2<-match(max(score2),score2)
    k=ifelse(score1[k1]<score2[k2], -1*(k2-1), (k1-1))
    max=ifelse(score1[k1]<score2[k2], score2[k2], score1[k1])
    list(k=k,max=max)
}