\name{plotMotifLogoStackWithTree}
\alias{plotMotifLogoStackWithTree}
\title{plot sequence logos stack with hierarchical cluster tree}
\description{
	plot sequence logos stack with hierarchical cluster tree
}
\usage{
plotMotifLogoStackWithTree(pfms, hc, treewidth=1/8, trueDist=FALSE, ...)
}
\arguments{
	\item{pfms}{a list of position frequency matrices, pfms must be a list of class pfm}
	\item{hc}{an object of the type produced by stats::hclust}
	\item{treewidth}{the width to show tree}
	\item{trueDist}{logical flags to use hclust height or not.}
	\item{...}{other parameters can be passed to plotMotifLogo function}
}
\value{
	none
}
\examples{
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
	pcm3<-matrix(c(50,50,0,50,
				   0,100,0,0,
				   0,50,50,0,
				   0,0,0,100,
				   50,50,0,0,
				   0,0,50,50), nrow=4)
	rownames(pcm1)<-c("A","C","G","T")
	rownames(pcm2)<-c("A","C","G","T")
	rownames(pcm3)<-c("A","C","G","T")
	pfms<-list(p1=new("pfm",mat=pcm2pfm(pcm1),name="m1"),
			   p2=new("pfm",mat=pcm2pfm(pcm2),name="m2"),
			   p3=new("pfm",mat=pcm2pfm(pcm3),name="m3"))
	pfms<-DNAmotifAlignment(pfms)
	dis<-function(mj,nj){
		sum((mj-nj)^2)
	}
	mdis<-function(m,n){
		d<-0
		for(j in 1:ncol(m)) d<-d+dis(m[,j],n[,j])
		d
	}
	m<-matrix(nrow=3,ncol=3,dimnames=list(c("p1","p2","p3"),c("p1","p2","p3")))
	for(i in 1:3)
		for(j in 1:3)
			m[i,j]<-mdis(pfms[[i]]@mat,pfms[[j]]@mat)
	d<-dist(m)
	hc<-hclust(d)
	plotMotifLogoStackWithTree(pfms, hc=hc)
}