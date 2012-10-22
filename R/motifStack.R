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
        }
		if(xaxis) plotXaxis(pfm, p)
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

DNAmotifAlignment<-function(pfms, threshold=0.4, minimalConsensus=0, rcpostfix="(RC)", revcomp=rep(TRUE, length(pfms))){
    if(length(pfms)<2) stop("less than 2 motifs")
    if(length(revcomp)!=length(pfms)) stop("length of revcomp and pfms is not identical")
    lapply(pfms,function(.ele){
           if(class(.ele)!="pfm") stop("pfms must be a list of class pfm")
           if(.ele@alphabet!="DNA") stop("the alphabet of pfm must be DNA")
           })
    pfmcopy<-list(length=length(pfms))
    pfmcopy[[1]]<-pfms[[1]]
    for(i in 2:length(pfms)){
        pfmcopy[[i]]<-pfms[[i]]
        pfmcopy<-UngappedAlignment(pfmcopy, i, threshold, minimalConsensus, rcpostfix, revcomp[i])
    }
    l<-unlist(lapply(pfmcopy,function(.ele) ncol(.ele@mat)))
    l<-max(l)-l
    for(i in 1:length(pfmcopy)){
        pfmcopy[[i]]<-addBlank(pfmcopy[[i]],l[i],TRUE)
    }
    pfmcopy
}

###############################################################################
######## plot motif logo without plot.new
######## to be used to create a better view of stack, eg. radial sty,
###############################################################################
plotMotifLogoA<-function(pfm, font="Helvetica-Bold"){
	if (class(pfm) != "pfm"){
		stop("pfms must be a list of class pfm")
	}
	rname<-rownames(pfm@mat)
	if(is.null(rname))
    stop("pfm rowname is empty")
# get ghostscript path
	gscmd <- Sys.getenv("R_GSCMD")
	npos<-ncol(pfm@mat)
	ncha<-nrow(pfm@mat) 
	symbols<-list()
	for(i in 1:ncha){
		ps<-paste("%!PS\n/",font," findfont\n72 scalefont\n",
				  motifStack:::hex2psrgb(pfm@color[i])," setrgbcolor\nsetfont\nnewpath\n0 0 moveto\n(",
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
	ic<-motifStack:::getIC(pfm)
	ie<-motifStack:::getIE(pfm@mat)
	dw<-1/npos
	x.pos<-0
	for(j in 1:npos){
		column<-pfm@mat[,j]
		heights<-column*ic[j]/ie
		id<-order(heights)
		
		y.pos<-0
		for(i in 1:ncha){
			h<-heights[id[i]]
			if(h>0) grid.draw(grImport::pictureGrob(symbols[[id[i]]],x.pos,y.pos,dw,h,just=c(0,0),distort=TRUE))
			y.pos<-y.pos+h
		}
		x.pos<-x.pos+dw
	}
}

###############################################################################
######## radial style stack
######## 
###############################################################################
plotMotifStackWithRadialPhylog <- function (phylog, pfms=NULL,  
circle=1, cleaves=1, cnodes=0, 
labels.leaves=names(phylog$leaves), clabel.leaves=1,
labels.nodes=names(phylog$nodes), clabel.nodes=0, 
draw.box=FALSE, 
col.leaves=rep("black", length(labels.leaves)), 
col.leaves.bg=NULL, col.leaves.bg.alpha=1,
col.bg=NULL, col.bg.alpha=1, 
col.inner.label.circle=NULL, col.outer.label.circle=NULL, outer.label.circle.width="default",
clockwise =FALSE, init.angle=if(clockwise) 90 else 0, 
angle=360) 
{
	if (!inherits(phylog, "phylog")) 
    stop("Non convenient data")
	leaves.number <- length(phylog$leaves)
	checkLength <- function(tobechecked){
		!((length(tobechecked)>=leaves.number)|is.null(tobechecked))
	}
	for(tobechecked in c("col.leaves", "col.leaves.bg", "col.bg", "col.inner.label.circle", "col.outer.label.circle")){
		if(checkLength(eval(as.symbol(tobechecked)))) stop(paste("the length of", tobechecked, "should be same as the length of leaves"))
	}
	leaves.names <- names(phylog$leaves)
	nodes.number <- length(phylog$nodes)
	nodes.names <- names(phylog$nodes)
	if (length(labels.leaves) != leaves.number) 
    labels.leaves <- names(phylog$leaves)
	if (length(labels.nodes) != nodes.number) 
    labels.nodes <- names(phylog$nodes)
	if (circle < 0) 
    stop("'circle': non convenient value")
	leaves.car <- gsub("[_]", " ", labels.leaves)
	nodes.car <- gsub("[_]", " ", labels.nodes)
	opar <- par(mar = par("mar"), srt = par("srt"))
	on.exit(par(opar))
	par(mar = c(0.1, 0.1, 0.1, 0.1))
	dis <- phylog$droot
	dis <- dis/max(dis)
	rayon <- circle
	dis <- dis * rayon
	dist.leaves <- dis[leaves.names]
	dist.nodes <- dis[nodes.names]
	asp <- c(1, 1)
	if(is.null(pfms)){
		plot.default(0, 0, type = "n", asp = 1, xlab = "", ylab = "", 
					 xaxt = "n", yaxt = "n", xlim = c(-2, 2), ylim = c(-2, 2), 
					 xaxs = "i", yaxs = "i", frame.plot = FALSE)
	}else{
		pin <- par("pin")
		if (pin[1L] > pin[2L]) asp <- c(pin[2L]/pin[1L], 1)
		else asp <- c(1, pin[1L]/pin[2L])
		plot.default(0, 0, type = "n", asp=1, xlab = "", ylab = "", 
					 xaxt = "n", yaxt = "n", xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5), 
					 xaxs = "i", yaxs = "i", frame.plot = FALSE)
	}
	d.rayon <- rayon/(nodes.number - 1)
	if(outer.label.circle.width=="default") outer.label.circle.width<-d.rayon
	twopi <- if (clockwise)
    -2 * pi
	else 2 * pi
	alpha <- twopi * angle * (1:leaves.number)/leaves.number/360 + init.angle * pi/180
	names(alpha) <- leaves.names
	x <- dist.leaves * cos(alpha)
	y <- dist.leaves * sin(alpha)
	xcar <- (rayon + d.rayon) * cos(alpha)
	ycar <- (rayon + d.rayon) * sin(alpha)
##for logos position
	if(!is.null(pfms)){
		beta <- alpha * 180 / pi 
		vpheight <- max(unlist(lapply(leaves.names, strheight, units="figure")))
		vpheight <- vpheight * asp[2L]
		mwidth <- max(unlist(lapply(pfms, function(.ele) ncol(.ele@mat))))
		vpwidth <- vpheight * mwidth / 2
		rayonWidth <- max(unlist(lapply(leaves.names, strwidth, units="user"))) * par("cex") * clabel.leaves
		xm <- (rayon + d.rayon + rayonWidth + 2.5*vpwidth) * cos(alpha) * asp[1L] / 5 + 0.5
		ym <- (rayon + d.rayon + rayonWidth + 2.5*vpwidth) * sin(alpha) * asp[2L] / 5 + 0.5
		if((max(xm)>1-vpwidth) | (min(xm)<vpwidth) | (max(ym)>1-vpheight) | (min(ym)<vpheight) ){
			if(interactive()){
				msg <- "Sequence logo will be drawn out of canvas. continue? (Y/n): "
				repeat{
					cat(msg)
					answer <- readLines(n=1)
					if(answer %in% c("Y","y","N","n"))
					break
					if(answer == ""){
						answer <- "y"
						break
					}  
				}
				tolower(answer)
				if(answer == "n"){
					return()
				}
			}
		}
	}
##for plot background
	if(!is.null(col.bg)) col.bg <- highlightCol(col.bg, col.bg.alpha)
	if(!is.null(col.leaves.bg)) col.leaves.bg <- highlightCol(col.leaves.bg, col.leaves.bg.alpha)
	gamma <- twopi * angle * ((1:(leaves.number+1))-0.5)/leaves.number/360 + init.angle * pi/180
	n <- max(2, floor(200*360/leaves.number))
	plotBgArc <- function(i, r, bgcol) {
		t2xy <- function(t) list(x=r*cos(t), y=r*sin(t))
		P <- t2xy(seq.int(gamma[i], gamma[i+1], length.out=n))
		polygon(c(P$x, 0), c(P$y, 0), border=bgcol, col=bgcol)
	}
	plotBgArc <- function(r,bgcol,inr){
		t2xy <- function(rx,t) list(x=rx*cos(t), y=rx*sin(t))
		oldcol <- bgcol[1]
		start <- 1
		icnt <- 1
		bgcol[length(bgcol)+1]<-NA
		for(i in 1:leaves.number){
			oldcol <- bgcol[i]
			if(bgcol[i+1]!=oldcol){
				P <- t2xy(r, seq.int(gamma[start], gamma[start+icnt], length.out=n*icnt))
				polygon(c(P$x, 0), c(P$y, 0), border=bgcol[i], col=bgcol[i])
				start <- i+1
				icnt <- 1
			} else {
				icnt <- icnt + 1
			}
		}
		if(inr!=0){
			P <- t2xy(inr, seq.int(0, twopi, length.out=n*leaves.number))
			polygon(c(P$x, 0), c(P$y, 0), border="white", col="white")
		}
	}
	if (clabel.leaves > 0) {
		if(!is.null(col.outer.label.circle)) ##plot outer.label.circle
		plotBgArc(rayon+d.rayon+rayonWidth+outer.label.circle.width, col.outer.label.circle, rayon+d.rayon+rayonWidth)
		if(!is.null(col.leaves.bg)) ##plot leaves bg
		plotBgArc(rayon+d.rayon+rayonWidth, col.leaves.bg, rayon+d.rayon)
		if(!is.null(col.inner.label.circle)) #plot inner.label.circle
		plotBgArc(rayon+d.rayon, col.inner.label.circle, rayon) 
		if(!is.null(col.bg)) ##plot center bg
		plotBgArc(rayon, col.bg, 0)
		for (i in 1:leaves.number) {
			segments(xcar[i], ycar[i], x[i], y[i], col = grey(0.7))
		}
		tmpout <- lapply(1:leaves.number,function(i) {
						 par(srt = alpha[i] * 180/pi)
						 text(xcar[i], ycar[i], leaves.car[i], adj = 0, col=col.leaves[i], cex = par("cex") * 
							  clabel.leaves)
						 segments(xcar[i], ycar[i], x[i], y[i], col = grey(0.7))
						 if(!is.null(pfms)){
						 pushViewport(viewport(x=xm[i], y=ym[i], width=vpwidth, height=vpheight, angle=beta[i]))
						 if(!is.null(pfms[[i]])){
						 plotMotifLogoA(pfms[[i]])
						 }
						 popViewport()
						 }
						 })
	}
	if (cleaves > 0) {
		for (i in 1:leaves.number) points(x[i], y[i], pch = 21, col=col.leaves[i],
										  bg = col.leaves[i], cex = par("cex") * cleaves)
	}
	ang <- rep(0, length(dist.nodes))
	names(ang) <- names(dist.nodes)
	ang <- c(alpha, ang)
	for (i in 1:length(phylog$parts)) {
		w <- phylog$parts[[i]]
		but <- names(phylog$parts)[i]
		ang[but] <- mean(ang[w])
		b <- range(ang[w])
		a.seq <- c(seq(b[1], b[2], by = pi/180), b[2])
		lines(dis[but] * cos(a.seq), dis[but] * sin(a.seq), col="#222222")
		x1 <- dis[w] * cos(ang[w])
		y1 <- dis[w] * sin(ang[w])
		x2 <- dis[but] * cos(ang[w])
		y2 <- dis[but] * sin(ang[w])
		segments(x1, y1, x2, y2, col="#222222")
	}
	if (cnodes > 0) {
		for (i in 1:length(phylog$parts)) {
			w <- phylog$parts[[i]]
			but <- names(phylog$parts)[i]
			ang[but] <- mean(ang[w])
			points(dis[but] * cos(ang[but]), dis[but] * sin(ang[but]), 
				   pch = 21, bg = "white", cex = par("cex") * cnodes)
		}
	}
	points(0, 0, pch = 21, cex = par("cex") * 2, bg = "red")
	if (clabel.nodes > 0) {
		delta <- strwidth(as.character(length(dist.nodes)), cex = par("cex") * 
						  clabel.nodes)
		for (j in 1:length(dist.nodes)) {
			i <- names(dist.nodes)[j]
			par(srt = (ang[i] * 360/2/pi + 90))
			x1 <- dis[i] * cos(ang[i])
			y1 <- dis[i] * sin(ang[i])
			symbols(x1, y1, delta, bg = "white", add = TRUE, 
					inches = FALSE)
			text(x1, y1, nodes.car[j], adj = 0.5, cex = par("cex") * 
				 clabel.nodes)
		}
	}
	if (draw.box) 
    box()
	return(invisible())
}