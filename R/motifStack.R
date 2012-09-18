plotMotifLogo<-function(pfm, motifName, p=rep(0.25, 4), font="Helvetica-Bold", 
colset=c("#00811B","#2000C7","#FFB32C","#D00001"), 
xaxis=TRUE,yaxis=TRUE,xlab="position",ylab="bits",
xlcex=1.2, ylcex=1.2, ncex=1.2){
    if (class(pfm) == "data.frame"){
        pfm <- as.matrix(pfm)
    }else if (class(pfm) != "matrix"){
        stop("pfm must be of class matrix or data.frame")
    }
    if (any(abs(1 - apply(pfm,2,sum)) > 0.01))
    stop("Columns of pfm must add up to 1.0")
    if(class(colset)!="character")
    stop("colset must be of class character")
    if (length(colset)!=nrow(pfm))
    stop(paste("colset length and pfm row numbers different",length(colset),"!=",nrow(pfm)))
    rname<-rownames(pfm)
    if(is.null(rname))
    stop("pfm rowname is empty")
# get ghostscript path
    gscmd <- Sys.getenv("R_GSCMD")
    if(is.null(gscmd) || !nzchar(gscmd)) {
        gscmd <- switch(.Platform$OS.type,
                        unix = "gs",
                        windows = {
                        poss <- Sys.which(c("gswin64c.exe",
                                            "gswin32c.exe"))
                        poss <- poss[nzchar(poss)]
                        if (length(poss)) gscmd <- poss
                        })
    }
    if(is.null(gscmd) || !nzchar(gscmd)){
		plot.new()
		text(.5,.5,"Can not locate ghostscript.\nFor windows, please install ghostscript first and then try:\nSys.setenv(R_GSCMD=\"\\\"C:\\\\Path\\\\to\\\\gsbin\\\\gswin32c.exe\\\"\")")
    }else{
        npos<-ncol(pfm)
        ncha<-nrow(pfm) 
        symbols<-list()
        for(i in 1:ncha){
            ps<-paste("%!PS\n/",font," findfont\n72 scalefont\n",
                      hex2psrgb(colset[i])," setrgbcolor\nsetfont\nnewpath\n0 0 moveto\n(",
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
    #calculate postion of each symbol and plot
        plot.new()
        
        ic<-getIC(pfm, p)
        ie<-getIE(pfm)
        dw<-1/npos
        x.pos<-0
        for(j in 1:npos){
            column<-pfm[,j]
            heights<-column*ic[j]/ie
            id<-order(heights)
            
            y.pos<-0
            for(i in 1:ncha){
                h<-heights[id[i]]
                if(h>0) grImport::picture(symbols[[id[i]]],x.pos,y.pos,x.pos+dw,y.pos+h)
                y.pos<-y.pos+h
            }
            x.pos<-x.pos+dw
            if(xaxis) plotXaxis(pfm, p)
        }
        if(yaxis) plotYaxis(pfm)
        if(!is.na(xlab)) mtext(xlab,1,line=2,cex=xlcex)
        if(!is.na(ylab)) mtext(ylab,2,line=2,cex=ylcex)
        if(!missing(motifName)) mtext(motifName,3,line=0,cex=ncex)
    }
}

plotXaxis<-function(pfm, p=rep(0.25, 4)){
    npos<-ncol(pfm)
    ic<-getIC(pfm,p)
    dw<-1/npos
    at<-c()
    label<-c()
    j=1
    for(i in 1:npos){
        if(ic[i]>0){
            at<-c(at,dw*(i-0.5))
            label<-c(label,j)
            j<-j+1
        }
    }
    axis(1,at=at,labels=label)
}

plotYaxis<-function(pfm){
    ie<-getIE(pfm)
    majorat<-seq(0,floor(ie)/ceiling(ie),by=1/ceiling(ie))
    majorlab<-0:floor(ie)
    axis(2,at=majorat,labels=majorlab)
    minorat<-seq(0,ie,by=1/(ceiling(ie)*5))
    minorat<-minorat[!(minorat %in% majorat)]
    axis(2,at=minorat,tcl=par("tcl")*0.5,labels=FALSE)
}

plotMotifLogoStack<-function(pfms, ...){
    n<-length(pfms)
    lapply(pfms,function(.ele){
           if(class(.ele)!="pfm") stop("pfms must be a list of class pfm")
           })
    opar<-par(mfrow=c(n,1),mar=c(3.5,3.5,1.5,0.5))
    for(i in 1:(n-1)){
        plot(pfms[[n-i+1]],xlab=NA, ...)
    }
    plot(pfms[[1]], ...)
    par(opar)
}

plotMotifLogoStackWithTree<-function(pfms, hc, treewidth=1/8, trueDist=FALSE, ...){
    n<-length(pfms)
    lapply(pfms,function(.ele){
           if(class(.ele)!="pfm") stop("pfms must be a list of class pfm")
           })
    if(class(hc)!="hclust") stop("hc class must be hclust")
    if(treewidth>0.5) stop("treewidth can not greater than 0.5")
    layout(matrix(c(rep(1,n),rep(2:(n+1),ceiling(1/treewidth)-1)),nrow=n,ncol=ceiling(1/treewidth)))
#plot tree
    opar<-par(mar=c(0,0,0,0))
    plot.new()
    if(trueDist) h <- hc$height / max(hc$height) / 1.05
    else h<- seq(0.01, 0.95, length.out=length(hc$height))[order(hc$height)]
    m <- hc$merge
    o <- hc$order
    
    m[m > 0] <- n + m[m > 0] 
    m[m < 0] <- abs(m[m < 0])
    
    dist <- matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c("x", "y"))) 
    dist[1:n, 2] <- 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)
    
    for(i in 1:nrow(m)){
        dist[n + i, 2] <- (dist[m[i, 1], 2] + dist[m[i, 2], 2]) / 2
        dist[n + i, 1] <- h[i]
    }
    
    trx <- function(x) 0.9*treewidth*(1 - x)+ 0.01
    for(i in 1:nrow(dist)){
        dist[i, 1] <- trx(dist[i, 1])
        dist[i, 2] <- dist[i, 2] + 1/(4*n)
    }
    
    draw_connection <- function(x1, x2, y1, y2, x){
        grid.lines(x = c(x1, x), y = c(y1, y1))
        grid.lines(x = c(x2, x), y = c(y2, y2))
        grid.lines(x = c(x, x), y = c(y1, y2))
    }
    
    for(i in 1:nrow(m)){
        draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], trx(h[i]))
    }
    
#plot logo
    par(mar=c(3.5,3.5,1.5,0.5))
    for(i in 1:(n-1)){
        plot(pfms[[n-i+1]],xlab=NA, ...)
    }
    plot(pfms[[1]], ...)
    par(opar)
}

DNAmotifAlignment<-function(pfms, threshold=0.4){
    if(length(pfms)<2) stop("less than 2 motifs")
    lapply(pfms,function(.ele){
           if(class(.ele)!="pfm") stop("pfms must be a list of class pfm")
           if(.ele@alphabet!="DNA") stop("the alphabet of pfm must be DNA")
           })
    pfmcopy<-list(length=length(pfms))
    pfmcopy[[1]]<-pfms[[1]]
    for(i in 2:length(pfms)){
        pfmcopy[[i]]<-pfms[[i]]
        pfmcopy<-UngappedAlignment(pfmcopy, i, threshold)
    }
    l<-unlist(lapply(pfmcopy,function(.ele) ncol(.ele@mat)))
    l<-max(l)-l
    for(i in 1:length(pfmcopy)){
        pfmcopy[[i]]<-addBlank(pfmcopy[[i]],l[i],TRUE)
    }
    pfmcopy
}