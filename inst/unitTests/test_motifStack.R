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
	pfm.gat1.jaspar = new("pfm", mat=matrix(c(0.1287129, 0.4356436, 0.1881188, 0.2475248, 
	                                          0.17, 0.42, 0.05, 0.36, 
	                                          0, 0, 1, 0, 
	                                          1, 0, 0, 0, 
	                                          0, 0, 0, 1, 
	                                          1, 0, 0, 0, 
	                                          1, 0, 0, 0, 
	                                          0.03030303, 0.09090909, 0.8484848, 0.03030303), 
	                                        nrow=4, 
	                                        dimnames=list(c("A","C","G","T"))), 
	                      name="GAT1-JASPAR")
	pfm.gat1.scertf = new("pfm", mat=matrix(c(0.01010101, 0.4646465, 0.2424242, 0.2828283, 
	                                          0.01, 0.85, 0.13, 0.01, 
	                                          0.01, 0.01, 0.01, 0.97, 
	                                          0.03, 0.01, 0.01, 0.95, 
	                                          0.97, 0.01, 0.01, 0.01, 
	                                          0.01, 0.01, 0.01, 0.97, 
	                                          0.01, 0.97, 0.01, 0.01, 
	                                          0.2244898, 0.01020408, 0.6632653, 0.1020408, 
	                                          0.02040816, 0.3163265, 0.4285714, 0.2346939), 
	                                        nrow=4, 
	                                        dimnames=list(c("A","C","G","T"))), 
	                      name="GAT1-ScerTF")
	pfm.gat1.uniprobe = new("pfm", mat=matrix(c(0.1844893, 0.3204164, 0.2288104, 0.2662839, 
	                                            0.1939449, 0.2938899, 0.2340073, 0.2781579, 
	                                            0.3622495, 0.168747, 0.2328523, 0.2361512, 
	                                            0.25089, 0.1950942, 0.1270668, 0.426949, 
	                                            0.1871758, 0.3151092, 0.1139613, 0.3837536, 
	                                            0.1168009, 0.2115293, 0.188105, 0.4835648, 
	                                            0.007896063, 0.4572002, 0.2442517, 0.290652, 
	                                            0.001486376, 0.8602653, 0.1327394, 0.005508941, 
	                                            0.003332312, 0.002954489, 0.002197378, 0.9915158, 
	                                            0.03030772, 0.0005158365, 0.0008668447, 0.9683096, 
	                                            0.9921533, 0.001267418, 0.004773612, 0.001805671, 
	                                            0.001866145, 0.002382037, 0.002223528, 0.9935283, 
	                                            0.00265525, 0.9924961, 0.001731091, 0.00311756, 
	                                            0.225576, 0.008512981, 0.6597229, 0.1061882, 
	                                            0.02166241, 0.3155481, 0.4275484, 0.2352412, 
	                                            0.2674831, 0.1705538, 0.4457152, 0.116248, 
	                                            0.2924191, 0.1403103, 0.2492224, 0.3180481, 
	                                            0.2725799, 0.1674882, 0.2519977, 0.3079342, 
	                                            0.25228, 0.2583783, 0.1255863, 0.3637555, 
	                                            0.2920132, 0.2484219, 0.1864655, 0.2730994), 
	                                          nrow=4, 
	                                          dimnames=list(c("A","C","G","T"))), 
	                        name="GAT1-UniPROBE")
	pfms.revcomp.test<-DNAmotifAlignment(c(pfm.gat1.jaspar, pfm.gat1.scertf, pfm.gat1.uniprobe), 
							 revcomp=c(TRUE, TRUE, FALSE))
	checkEquals("GAT1-ScerTF(RC)", pfms.revcomp.test[[2]]@name)
	checkEquals("GAT1-UniPROBE", pfms.revcomp.test[[3]]@name)
}

test_hex2psrgb<-function(){
	checkEquals("1 0 0", motifStack:::hex2psrgb("red"))
	checkEquals("0 1 0", motifStack:::hex2psrgb("#00FF00"))
	checkEquals("0 0 1", motifStack:::hex2psrgb("#0000FF"))
	checkEquals("0.8 0.8 0.8", motifStack:::hex2psrgb("#CCCCCC"))
}

test_addPseudolog2<-function(){
	checkEquals(.Machine$double.min.exp, motifStack:::addPseudolog2(0))
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