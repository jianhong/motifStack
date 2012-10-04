test_DNAmotifAlignment<-function(){
	pcm1<-matrix(c(0,50,0,50,
				   100,0,0,0,
				   0,100,0,0,
				   0,0,100,0,
				   0,0,0,100,
				   50,50,0,0,
				   0,0,50,50), nrow=4)
	pcm2<-matrix(c(50,50,0,0,
				   0,100,0,0,
				   0,50,50,0,
				   0,0,0,100,
				   50,50,0,0,
				   0,0,50,50), nrow=4)
	rownames(pcm1)<-c("A","C","G","T")
	rownames(pcm2)<-c("A","C","G","T")
	pfms<-list(p1=new("pfm",mat=pcm2pfm(pcm1),name="m1"),
			   p2=new("pfm",mat=pcm2pfm(pcm2),name="m2"))
	pfms<-DNAmotifAlignment(pfms)
	for(i in 1:length(pfms)){
		checkEquals(ncol(pfms[[i]]@mat),7)
	}
	blank=rep(.25,4)
	names(blank)=c("A","C","G","T")
	checkEquals(blank,pfms[[2]]@mat[,1])
}

test_hex2psrgb<-function(){
	checkEquals("1 0 0", motifStack:::hex2psrgb("red"))
	checkEquals("0 1 0", motifStack:::hex2psrgb("#00FF00"))
	checkEquals("0 0 1", motifStack:::hex2psrgb("#0000FF"))
	checkEquals("0.8 0.8 0.8", motifStack:::hex2psrgb("#CCCCCC"))
}

test_addPseudolog2<-function(){
	checkEquals(-10, motifStack:::addPseudolog2(0))
	checkEquals(0, motifStack:::addPseudolog2(1))
	checkEquals(1, motifStack:::addPseudolog2(2))
}

test_getIE<-function(){
	x<-matrix(1:8,nrow=4)
	checkEquals(2, motifStack:::getIE(x))
	x<-matrix(1:20,nrow=20)
	checkEquals(4.322, motifStack:::getIE(x), tolerance=0.01)
}

test_getICbyBase<-function(){
	checkEquals(c(0,0,0,2,2), motifStack:::getICbyBase(p=rep(0.25,4),c(0,0,0,1)))
	checkEquals(c(0,0,0.5,0.5,1.0), motifStack:::getICbyBase(p=rep(0.25,4),c(0,0,.5,.5)))
	checkEquals(c(0,0,0,0,0), motifStack:::getICbyBase(p=rep(0.25,4),c(.25,.25,.25,.25)))
	checkEquals(c(0,0,0,0,0), motifStack:::getICbyBase(p=c(0,0,0,1),c(0,0,0,1)))
}

test_getoffsetPosByIC<-function(){
	pcm1<-matrix(c(0,50,0,50,
				   100,0,0,0,
				   0,100,0,0,
				   0,0,100,0,
				   0,0,0,100,
				   50,50,0,0,
				   0,0,50,50), nrow=4)
	pcm2<-matrix(c(50,50,0,0,
				   0,100,0,0,
				   0,50,50,0,
				   0,0,0,100,
				   50,50,0,0,
				   0,0,50,50), nrow=4)
	rownames(pcm1)<-c("A","C","G","T")
	rownames(pcm2)<-c("A","C","G","T")
	pfms<-list(p1=new("pfm",mat=pcm2pfm(pcm1),name="m1"),
			   p2=new("pfm",mat=pcm2pfm(pcm2),name="m2"))
	offset<-motifStack:::getoffsetPosByIC(pfms[[1]],pfms[[2]],0.4)
	checkEquals(offset$k, 1)
	checkEquals(offset$max, 6)
}

test_getAlignedICWithoutGap<-function(){
	pcm1<-matrix(c(0,50,0,50,
				   100,0,0,0,
				   0,100,0,0,
				   0,0,100,0,
				   0,0,0,100,
				   50,50,0,0,
				   0,0,50,50), nrow=4)
	pcm2<-matrix(c(50,50,0,0,
				   0,100,0,0,
				   0,50,50,0,
				   0,0,0,100,
				   50,50,0,0,
				   0,0,50,50), nrow=4)
	rownames(pcm1)<-c("A","C","G","T")
	rownames(pcm2)<-c("A","C","G","T")
	pfms<-list(p1=new("pfm",mat=pcm2pfm(pcm1),name="m1"),
			   p2=new("pfm",mat=pcm2pfm(pcm2),name="m2"))
	offset<-motifStack:::getAlignedICWithoutGap(pfms[[1]],pfms[[2]],0.4)
	checkEquals(offset$offset, 1)
	checkEquals(offset$rev, FALSE)
}

test_getIC<-function(){
	pcm1<-matrix(c(0,50,0,50,
				   100,0,0,0,
				   0,100,0,0,
				   0,0,100,0,
				   0,0,0,100,
				   50,50,0,0,
				   0,0,50,50), nrow=4)
	checkEquals(c(1, 2, 2, 2, 2, 1, 1), getIC(pcm2pfm(pcm1),p=rep(0.25,4)))
}

test_matrixReverseComplement<-function(){
	pcm2<-matrix(c(50,50,0,0,
				   0,100,0,0,
				   0,50,50,0,
				   0,0,0,100,
				   50,50,0,0,
				   0,0,50,50), nrow=4)
	rownames(pcm2)<-c("A","C","G","T")
	p2<-new("pfm",mat=pcm2pfm(pcm2),name="m2")
	pr<-matrixReverseComplement(p2)
	checkEqualsNumeric(1, pr@mat[3,5])
}

test_plotMotifStackWithRadialPhylog<-function(){
}