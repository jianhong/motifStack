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
		key<-paste("x", ncha, font, paste(colset, collapse=""), paste(rname, collapse=""), sep="_")
		symbolsCache <- if(exists("tmp_motifStack_symbolsCache", where=".GlobalEnv")) get("tmp_motifStack_symbolsCache", pos=".GlobalEnv") else list()
		if(!is.null(symbolsCache[[key]])) symbols<-symbolsCache[[key]]
		else {
			symbols<-motifStack:::coloredSymbols(ncha, font, colset, rname)
			symbolsCache[[key]]<-symbols
			assign("tmp_motifStack_symbolsCache", symbolsCache, pos=".GlobalEnv")
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
	assign("tmp_motifStack_symbolsCache", list(), pos=".GlobalEnv")
    for(i in 1:(n-1)){
        plot(pfms[[n-i+1]],xlab=NA, ...)
    }
    plot(pfms[[1]], ...)
	rm(list="tmp_motifStack_symbolsCache", pos=".GlobalEnv")
    par(opar)
}

plotMotifLogoStackWithTree<-function(pfms, hc, treewidth=1/8, trueDist=FALSE, ...){
    n<-length(pfms)
    lapply(pfms,function(.ele){
           if(class(.ele)!="pfm") stop("pfms must be a list of class pfm")
           })
    if(class(hc)!="hclust") stop("hc class must be hclust")
    if(treewidth>0.5) stop("treewidth can not greater than 0.5")
    opar<-par(mar=c(0,0,0,0), mfrow=par("mfrow"))
    layout(matrix(c(rep(1,n),rep(2:(n+1),ceiling(1/treewidth)-1)),nrow=n,ncol=ceiling(1/treewidth)))
#plot tree
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
	assign("tmp_motifStack_symbolsCache", list(), pos=".GlobalEnv")
    for(i in 1:(n-1)){
        plot(pfms[[n-i+1]],xlab=NA, ...)
    }
    plot(pfms[[1]], ...)
	rm(list="tmp_motifStack_symbolsCache", pos=".GlobalEnv")
    par(opar)
}

DNAmotifAlignment<-function(pfms, threshold=0.4, minimalConsensus=0, rcpostfix="(RC)", revcomp=rep(TRUE, length(pfms))){
    if(length(pfms)<2) stop("less than 2 motifs")
    if(length(revcomp)!=length(pfms)) stop("length of revcomp and pfms is not identical")
    lapply(pfms,function(.ele){
           if(class(.ele)!="pfm" && class(.ele)!="pcm") stop("pfms must be a list of class pfm")
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
	key<-paste("x", ncha, font, paste(pfm@color, collapse=""), paste(rname, collapse=""), sep="_")
	symbolsCache <- if(exists("tmp_motifStack_symbolsCache", where=".GlobalEnv")) get("tmp_motifStack_symbolsCache", pos=".GlobalEnv") else list()
	if(!is.null(symbolsCache[[key]])) symbols<-symbolsCache[[key]]
	else {
		symbols<-motifStack:::coloredSymbols(ncha, font, pfm@color, rname)
		symbolsCache[[key]]<-symbols
		assign("tmp_motifStack_symbolsCache", symbolsCache, pos=".GlobalEnv")
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
######## phylog style stack
######## 
###############################################################################

plotMotifStackWithPhylog <- function(phylog, pfms=NULL,
f.phylog = 0.3, f.logo = NULL, cleaves =1, cnodes =0,
labels.leaves = names(phylog$leaves), clabel.leaves=1,
labels.nodes = names(phylog$nodes), clabel.nodes = 0
){
	if(!inherits(phylog, "phylog")) stop("phylog must be an object of phylog")
	n<-length(pfms)
	lapply(pfms,function(.ele){
		   if(class(.ele)!="pfm") stop("pfms must be a list of class pfm")
		   })
	leaves.number <- length(phylog$leaves)
	leaves.names<- names(phylog$leaves)
	nodes.number <- length(phylog$nodes)
	nodes.names <- names(phylog$nodes)
	if (length(labels.leaves) != leaves.number) labels.leaves <- names(phylog$leaves)
	if (length(labels.nodes) != nodes.number) labels.nodes <- names(phylog$nodes)
	leaves.car <- gsub("[_]"," ",labels.leaves)
	nodes.car <- gsub("[_]"," ",labels.nodes)
	opar <- par(mar = c(0, 0, 0, 0))
	on.exit(par(opar))
	
	if (f.phylog < 0.05) f.phylog <- 0.05 
	if (f.phylog > 0.75) f.phylog <- 0.75
	
	maxx <- max(phylog$droot)
	plot.default(0, 0, type = "n", xlab = "", ylab = "", xaxt = "n", 
				 yaxt = "n", xlim = c(-maxx*0.15, maxx/f.phylog), ylim = c(-0, 1), xaxs = "i", 
				 yaxs = "i", frame.plot = FALSE)
	
	x.leaves <- phylog$droot[leaves.names]
	x.nodes <- phylog$droot[nodes.names]
	
	y <- (leaves.number:1)/(leaves.number + 1)
	names(y) <- leaves.names
	
	xcar <- maxx*1.05
	xx <- c(x.leaves, x.nodes)
	
	assign("tmp_motifStack_symbolsCache", list(), pos=".GlobalEnv")
	if(is.null(f.logo)){
		f.logo <- max(unlist(lapply(leaves.names, strwidth, units="figure", cex=clabel.leaves)))
		if(clabel.leaves>0) f.logo <- f.phylog+f.logo+.01
		else f.logo <- f.phylog+.05
	}else{
		if (f.logo > 0.85) f.logo <- 0.85
		if (f.logo < f.phylog) f.logo <- f.phylog+.05 
	}
	for (i in 1:leaves.number) {
		if(clabel.leaves>0) 
		text(xcar, y[i], leaves.car[i], adj = 0, cex = par("cex") * 
			 clabel.leaves)
		segments(xcar, y[i], xx[i], y[i], col = grey(0.7))
		if(!is.null(pfms)){
			vpheight <- strheight(leaves.names[i], units="figure")
			vpwidth <- vpheight * ncol(pfms[[i]]@mat) / 2
			pushViewport(viewport(x=f.logo, y=y[i], width=vpwidth, height=vpheight, just=c(0, .5)))
			if(!is.null(pfms[[i]])) plotMotifLogoA(pfms[[i]])
			popViewport()
		}
	}
	rm(list="tmp_motifStack_symbolsCache", pos=".GlobalEnv")
	
	yleaves <- y[1:leaves.number]
	xleaves <- xx[1:leaves.number]
	if (cleaves > 0) {
		for (i in 1:leaves.number) {
			points(xx[i], y[i], pch = 21, bg=1, cex = par("cex") * cleaves)
		}
	}
	yn <- rep(0, nodes.number)
	names(yn) <- nodes.names
	y <- c(y, yn)
	for (i in 1:length(phylog$parts)) {
		w <- phylog$parts[[i]]
		but <- names(phylog$parts)[i]
		y[but] <- mean(y[w])
		b <- range(y[w])
		segments(xx[but], b[1], xx[but], b[2])
		x1 <- xx[w]
		y1 <- y[w]
		x2 <- rep(xx[but], length(w))
		segments(x1, y1, x2, y1)
	}
	if (cnodes > 0) {
		for (i in nodes.names) {
			points(xx[i], y[i], pch = 21, bg="white", cex = cnodes)
		}
	}
	if (clabel.nodes > 0) {
		scatterutil.eti(xx[names(x.nodes)], y[names(x.nodes)], nodes.car, 
						clabel.nodes)
	}
	if (cleaves > 0) points(xleaves, yleaves, pch = 21, bg=1, cex = par("cex") * cleaves)
	return(invisible())
}

###############################################################################
######## radial style stack
######## 
###############################################################################
plotMotifStackWithRadialPhylog <- function (phylog, pfms=NULL,
circle=1, circle.motif=NA, cleaves=1, cnodes=0,
labels.leaves=names(phylog$leaves), clabel.leaves=1,
labels.nodes=names(phylog$nodes), clabel.nodes=0,
draw.box=FALSE,
col.leaves=rep("black", length(labels.leaves)),
col.leaves.bg=NULL, col.leaves.bg.alpha=1,
col.bg=NULL, col.bg.alpha=1,
col.inner.label.circle=NULL, inner.label.circle.width="default",
col.outer.label.circle=NULL, outer.label.circle.width="default",
clockwise =FALSE, init.angle=if(clockwise) 90 else 0,
angle=360, pfmNameSpliter=";", rcpostfix="(RC)", motifScale=c("linear","logarithmic"))
{
	if (!inherits(phylog, "phylog"))
    stop("Non convenient data")
	leaves.number <- length(phylog$leaves)
	checkLength <- function(tobechecked){
		!((length(tobechecked)>=leaves.number)||is.null(tobechecked))
	}
	checkNA <- function(tobechecked){
		if(is.null(tobechecked)) return(FALSE)
		return(any(is.na(tobechecked)))
	}
	for(tobechecked in c("col.leaves", "col.leaves.bg", "col.bg", "col.inner.label.circle", "col.outer.label.circle")){
		if(checkLength(eval(as.symbol(tobechecked)))) stop(paste("the length of", tobechecked, "should be same as the length of leaves"))
		if(checkNA(eval(as.symbol(tobechecked)))) stop(paste("contain NA in", tobechecked))
	}
    motifScale <- match.arg(motifScale)
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
	par(mar = c(0.1, 0.1, 0.1, 0.1), mfrow=c(1,1))
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
	    pin <- dev.size("in") #pin <- par("pin")
		if (pin[1L] > pin[2L]) asp <- c(pin[2L]/pin[1L], 1)
		else asp <- c(1, pin[1L]/pin[2L])
		plot.default(0, 0, type = "n", asp=1, xlab = "", ylab = "",
        xaxt = "n", yaxt = "n", xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5),
        xaxs = "i", yaxs = "i", frame.plot = FALSE)
	}
	d.rayon <- rayon/(nodes.number - 1)
	twopi <- if (clockwise)
    -2 * pi
	else 2 * pi
	alpha <- twopi * angle * (1:leaves.number)/leaves.number/360 + init.angle * pi/180
	names(alpha) <- leaves.names
	x <- dist.leaves * cos(alpha)
	y <- dist.leaves * sin(alpha)
	xcar <- (rayon + d.rayon) * cos(alpha)
	ycar <- (rayon + d.rayon) * sin(alpha)
	
	rayonWidth <- max(unlist(lapply(leaves.names, strwidth, units="user", cex=clabel.leaves)))
	circle.motif <- ifelse(is.na(circle.motif), rayon + d.rayon + rayonWidth, circle.motif)
    ##for logos position
    maxvpwidth <- 0
	if(!is.null(pfms)){
		beta <- alpha * 180 / pi
		vpheight <- strheight("ACGT", units="figure")
		vpheight <- vpheight * asp[2L]
		xm <- circle.motif * cos(alpha) * asp[1L] / 5 + 0.5
		ym <- circle.motif * sin(alpha) * asp[2L] / 5 + 0.5
		maxvpwidth <- max(unlist(lapply(pfms, function(.ele) ncol(.ele@mat)))) * vpheight * 2.5
	}
	if(inner.label.circle.width=="default") inner.label.circle.width <- rayonWidth/10
	if(outer.label.circle.width=="default") outer.label.circle.width <- rayonWidth/10
    ##for plot background
	if(!is.null(col.bg)) col.bg <- motifStack:::highlightCol(col.bg, col.bg.alpha)
	if(!is.null(col.leaves.bg)) col.leaves.bg <- motifStack:::highlightCol(col.leaves.bg, col.leaves.bg.alpha)
	gamma <- twopi * angle * ((1:(leaves.number+1))-0.5)/leaves.number/360 + init.angle * pi/180
	n <- max(2, floor(200*360/leaves.number))
	plotBgArc <- function(r,bgcol,inr){
		t2xy <- function(rx,t) list(x=rx*cos(t), y=rx*sin(t))
		oldcol <- bgcol[1]
		start <- 1
		icnt <- 1
		for(i in 1:leaves.number){
			oldcol <- bgcol[i]
			if(i==leaves.number || bgcol[i+1]!=oldcol){
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
		plotBgArc(circle.motif+maxvpwidth+outer.label.circle.width, col.outer.label.circle, circle.motif+maxvpwidth)
		if(!is.null(col.leaves.bg)) ##plot leaves bg
		plotBgArc(circle.motif, col.leaves.bg, rayon+d.rayon)
		if(!is.null(col.inner.label.circle)) #plot inner.label.circle
		plotBgArc(circle.motif, col.inner.label.circle, circle.motif - inner.label.circle.width)
		if(!is.null(col.bg)) ##plot center bg
		plotBgArc(mean(dist.leaves), col.bg, 0)
		for(i in 1:leaves.number) {
			par(srt = alpha[i] * 180/pi)
			text(xcar[i], ycar[i], leaves.car[i], adj = 0, col=col.leaves[i], cex = par("cex") *
            clabel.leaves)
			segments(xcar[i], ycar[i], x[i], y[i], col = grey(0.7))
		}
		
		assign("tmp_motifStack_symbolsCache", list(), pos=".GlobalEnv")
		if(!is.null(pfms)){
			##extract names
            for(metaChar in c("\\","$","*","+",".","?","[","]","^","{","}","|","(",")"))
            {
                rcpostfix <- gsub(metaChar,paste("\\",metaChar,sep=""),rcpostfix,fixed=TRUE)
            }
			pfmNames <- lapply(pfms, function(.ele) .ele@name)
			for(i in 1:length(pfmNames)){
				pfmname <- unlist(strsplit(pfmNames[[i]], pfmNameSpliter))
                pfmname <- gsub(paste(rcpostfix,"$",sep=""),"",pfmname)
				pfmIdx <- which(makeLeaveNames(labels.leaves) %in% makeLeaveNames(pfmname))
                if(length(pfmIdx)==0) 
                    pfmIdx <- which(makeLeaveNames(names(phylog$leaves)) 
                                    %in% makeLeaveNames(pfmname))
				if(length(pfmIdx)>0){
					vph <- ifelse(motifScale=="linear",
                                    vpheight*length(pfmname),
                                    vpheight*(1+log2(length(pfmname))))
					vpw <- vph * ncol(pfms[[i]]@mat) / 2
					vpd <- sqrt(vph*vph+vpw*vpw) / 2
					vpx <- median(xm[pfmIdx]) + vpd * cos(median(alpha[pfmIdx])) * asp[1L]
					vpy <- median(ym[pfmIdx]) + vpd * sin(median(alpha[pfmIdx])) * asp[2L]
					angle <- median(beta[pfmIdx])
					pushViewport(viewport(x=vpx, y=vpy, width=vpw, height=vph, angle=angle))
					plotMotifLogoA(pfms[[i]])
					popViewport()
				}else{
                    warning(paste("No leave named as ", paste(pfmname, collapse=", ")), sep="")
                }
			}
		}
		rm(list="tmp_motifStack_symbolsCache", pos=".GlobalEnv")
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


###############################################################################
####### motifstack
motifStack <-function(pfms, layout=c("stack", "treeview", "phylog", "radialPhylog"), ...){
	layout <- match.arg(layout, c("stack", "treeview", "phylog", "radialPhylog"), several.ok=TRUE)
	layout <- layout[1]
	if (any(unlist(lapply(pfms, function(.ele) !inherits(.ele, "pfm"))))) 
    stop("pfms must be a list of pfm objects")
	if (length(pfms)<2)
    stop("length of pfms less than 2")
	pfmList2matrixList <- function(pfms){
		m <- lapply(pfms, function(.ele) as(.ele, "matrix"))
		names(m) <- unlist(lapply(pfms, function(.ele) .ele@name))
		m
	}
	switch(layout,
		   stack = {
           plotMotifLogoStack(pfms, ...)
		   },
		   treeview = {
		   jaspar.scores <- MotIV::readDBScores(file.path(find.package("MotIV"), "extdata", "jaspar2010_PCC_SWU.scores"))
           d <- MotIV::motifDistances(pfmList2matrixList(pfms))
           hc <- MotIV::motifHclust(d)
           pfms <- pfms[hc$order]
           pfms <- DNAmotifAlignment(pfms)
           plotMotifLogoStackWithTree(pfms, hc=hc, ...)
		   },
		   phylog = {
		   jaspar.scores <- MotIV::readDBScores(file.path(find.package("MotIV"), "extdata", "jaspar2010_PCC_SWU.scores"))
           d <- MotIV::motifDistances(pfmList2matrixList(pfms))
           hc <- MotIV::motifHclust(d)
           pfms <- pfms[hc$order]
           pfms <- DNAmotifAlignment(pfms)
           phylog <- hclust2phylog(hc)
           plotMotifStackWithPhylog(phylog=phylog, pfms=pfms, ...)
		   },
		   radialPhylog = {
		   jaspar.scores <- MotIV::readDBScores(file.path(find.package("MotIV"), "extdata", "jaspar2010_PCC_SWU.scores"))
           d <- MotIV::motifDistances(pfmList2matrixList(pfms))
           hc <- MotIV::motifHclust(d)
           pfms <- pfms[hc$order]
           pfms <- DNAmotifAlignment(pfms)
           phylog <- hclust2phylog(hc)
           args <- list(phylog=phylog, pfms=pfms, ...)
           for(i in names(args)){
               if(i %in% c("col.leaves", "col.leaves.bg", "col.bg", "col.inner.label.circle", "col.outer.label.circle")){
                   args[[i]] <- args[[i]][hc$order]
               }
           }
           do.call(plotMotifStackWithRadialPhylog, args)
		   },
		   {
		   plotMotifLogoStack(pfms, ...)
		   }
		   )
}

###############################################################################
######## motif cloud
######## 
###############################################################################
motifSignature <- function(pfms, phylog, groupDistance, rcpostfix="(RC)", 
min.freq=2, trim=0.2, families=list()){
	if (!inherits(phylog, "phylog")) 
    stop("phylog should be an object of phylog of package ade4")
	leaves.number <- length(phylog$leaves)
	if (length(pfms)!=leaves.number)
    stop("length of pfms and leaves of phylog should be identical.")
	if (all(unlist(lapply(pfms, function(.ele) !inherits(.ele, "pfm"))))){
		pfms.class="pcm"
		if(any(unlist(lapply(pfms, function(.ele) !inherits(.ele, "pcm")))))
			stop("pfms should be a list of objects of pfm or pcm")
	}else{
		pfms.class="pfm"
		if(any(unlist(lapply(pfms, function(.ele) !inherits(.ele, "pfm")))))
			stop("pfms should be a list of objects of pfm or pcm")
	}
	if (missing(groupDistance))
		groupDistance <- max(phylog$droot)/10
#generate the new phylog from newick tree
	tree <- phylog$tre
	droot <- phylog$droot
	tree <- gsub(";$","",tree)
	str <- unlist(strsplit(gsub("^\\((.*?)\\)([^\\)\\(;,]+)$","\\1///\\2",tree),"///",fixed=T))[1]
	nodelist <- list()
	getNodelist <- function(str, nodelist=list()){
		ge <- gregexpr("\\(([^\\(\\),]+),([^\\(\\),]+)\\)([^\\(\\),]+)", str)
		if(ge[[1]][1]!=-1){
			start <- ge[[1]]
			stop <- attr(ge[[1]],"match.length")
			trec<-c()
			for(i in 1:length(start)){
				trec<- c(trec,substr(str, start[i],start[i]+stop[i]-1))
			}
			nodes <- do.call(rbind,strsplit(trec,"\\)|\\(|,"))[,2:4,drop=F]
			newnodes <- apply(nodes, 1, function(.ele) new("ouNode",left=.ele[1], right=.ele[2], parent=.ele[3], sizel=0, sizer=0))
			names(newnodes) <- nodes[,3]
			nodelist <- c(newnodes,nodelist)
			for(i in 1:length(trec)){
				str<-gsub(trec[i],gsub(".*?\\)(.*?)$","\\1",trec[i]),str,fixed=T)
			}
			Recall(str=str, nodelist=nodelist)
		}else{
			nodes <- unlist(strsplit(str, ","))
			Root <- new("ouNode", left=nodes[1], right=nodes[2], parent="Root", sizel=0, sizer=0)
			nodelist <- c(Root=Root, nodelist)
			return(nodelist)
		}
	}
	nodelist <- getNodelist(str)
	leaves <- names(phylog$leaves)
	buildTree <- function(nodelist, nodename, droot, dpar=0){
		nodelist[[nodename]]@distl <<- droot[nodelist[[nodename]]@left] - dpar
		nodelist[[nodename]]@distr <<- droot[nodelist[[nodename]]@right] - dpar
		if(nodelist[[nodename]]@left %in% leaves) nodelist[[nodename]]@sizel <<- 1
		if(nodelist[[nodename]]@right %in% leaves) nodelist[[nodename]]@sizer <<- 1
		if(!is.null(nodelist[[nodelist[[nodename]]@left]])){
			buildTree(nodelist, nodelist[[nodename]]@left, droot, droot[nodelist[[nodename]]@left])
		}
		if(!is.null(nodelist[[nodelist[[nodename]]@right]])){
			buildTree(nodelist, nodelist[[nodename]]@right, droot, droot[nodelist[[nodename]]@right])
		}
	}
	buildTree(nodelist, "Root", droot, dpar=0)
	
	mergeNodes <- function(nodelist, leaves, groupDistance){
		l <- length(nodelist)
		for(i in 1:l){
			nodename <- names(nodelist)[i]
			currNode <- nodelist[[nodename]]
			if(currNode@left %in% leaves && currNode@right %in% leaves){
				if(currNode@distl < groupDistance && currNode@distr < groupDistance){
					parNodeInfo <- getParentNode(nodelist, nodename)
					if(!is.null(parNodeInfo)[1]){
						if(parNodeInfo[2]=="left") {
							nodelist[[parNodeInfo[1]]]@sizel <- nodelist[[parNodeInfo[1]]]@sizel + max(currNode@sizel, 1) + max(currNode@sizer, 1)
							nodelist[[parNodeInfo[1]]]@distl <- nodelist[[parNodeInfo[1]]]@distl + (currNode@distl + currNode@distr)/2
							nodelist[[parNodeInfo[1]]]@left <- paste(currNode@left, currNode@right, sep=";")
							leaves <- c(leaves[!leaves %in% c(currNode@left, currNode@right)], paste(currNode@left, currNode@right, sep=";"))
							nodelist[[nodename]] <- new("ouNode",left="NULL", right="NULL", parent="NULL", sizel=0, sizer=0)
						}
						else if(parNodeInfo[2]=="right") {
							nodelist[[parNodeInfo[1]]]@sizer <- nodelist[[parNodeInfo[1]]]@sizer + max(currNode@sizel, 1) + max(currNode@sizer, 1)
							nodelist[[parNodeInfo[1]]]@distr <- nodelist[[parNodeInfo[1]]]@distr + (currNode@distl + currNode@distr)/2
							nodelist[[parNodeInfo[1]]]@right <- paste(currNode@left, currNode@right, sep=";")
							leaves <- c(leaves[!leaves %in% c(currNode@left, currNode@right)], paste(currNode@left, currNode@right, sep=";"))
							nodelist[[nodename]] <- new("ouNode",left="NULL", right="NULL", parent="NULL", sizel=0, sizer=0)
						}
					}
				}
			}
		}
		sel <- unlist(lapply(nodelist, function(.ele) .ele@parent=="NULL"))
		nodelist <- nodelist[!sel]
		nodelist <<- nodelist
		leaves <<- leaves
		if(length(nodelist)<l) mergeNodes(nodelist, leaves, groupDistance)
	}
	
	getSignature <- function(pfms, pfmnames, rcpostfix, pfms.class){
		.pfmnames <- unlist(strsplit(pfmnames, ";"))
		.pfms <- lapply(.pfmnames, function(.name, pfms) pfms[[getPFMid(pfms, .name, rcpostfix=rcpostfix)]], pfms)
		if(length(.pfms)>1){
			.pfms <- DNAmotifAlignment(.pfms, rcpostfix="")
			.rown <- rownames(.pfms[[1]]@mat)
			.mats <- lapply(.rown, function(i, .pfms){
								do.call(rbind,lapply(.pfms, function(.ele, i) .ele@mat[i, , drop=F], i))
							}, .pfms)
			if(pfms.class=="pfms"){
				.mats <- do.call(rbind, lapply(.mats, colMeans))
			}
			else{
				.mats <- do.call(rbind, lapply(.mats, colSums))
				.mats <- pcm2pfm(.mats)
			}
			rownames(.mats) <- .rown
			colnames(.mats) <- 1:ncol(.mats)
		}else{
			if(pfms.class=="pfms") .mats <- .pfms[[1]]@mat
			else .mats <- pcm2pfm(.pfms[[1]]@mat)
		}
		new("pfm", mat=.mats, name=pfmnames, alphabet=.pfms[[1]]@alphabet, color=.pfms[[1]]@color, background=.pfms[[1]]@background)
	}
	
	mergeNodes(nodelist, leaves, groupDistance)
#get signatures and frequences (count of frequences should greater than min.freq)
	filterFamilies <- function(pfmnames, size, families){
		if(length(families)>0){
			pfmnames <- unlist(strsplit(pfmnames, ";"))
			if(length(pfmnames)!=size) stop("pfm names contain ';'")
			for(j in 1:length(families)){
				inters <- intersect(families[[j]], pfmnames)
				if(length(inters)>1){
					pfmnames <- pfmnames[!pfmnames %in% inters[2:length(inters)]]
				}
			}
			size <- length(pfmnames)
			pfmnames <- paste(pfmnames, sep="", collapse=";")
		}
		list(name=pfmnames,size=size)
	}
	signatures <- list()
	freq <- c()
	for(i in 1:length(nodelist)){
		if(nodelist[[i]]@sizel>=min.freq){
			ftmp <- filterFamilies(nodelist[[i]]@left, nodelist[[i]]@sizel, families)
			nodelist[[i]]@left <- ftmp[["name"]]
			nodelist[[i]]@sizel <- ftmp[["size"]]
			if(ftmp[["size"]]>=min.freq){
				signatures <- c(signatures, getSignature(pfms, nodelist[[i]]@left, rcpostfix, pfms.class=pfms.class))
				freq <- c(freq, nodelist[[i]]@sizel)
			}
		}
		if(nodelist[[i]]@sizer>=min.freq){
			ftmp <- filterFamilies(nodelist[[i]]@right, nodelist[[i]]@sizer, families)
			nodelist[[i]]@right <- ftmp[["name"]]
			nodelist[[i]]@sizer <- ftmp[["size"]]
			if(ftmp[["size"]]>=min.freq){
				signatures <- c(signatures, getSignature(pfms, nodelist[[i]]@right, rcpostfix, pfms.class=pfms.class))
				freq <- c(freq, nodelist[[i]]@sizer)
			}
		}
	}
	
#trime signatures
	signatures <- lapply(signatures, trimMotif, trim)
	ord <- unlist(lapply(signatures, function(.ele){inherits(.ele,"pfm")}))
	signatures <- signatures[ord]
	freq <- freq[ord]
#sort signatures
	ord <- order(freq, decreasing=TRUE)
	signatures <- signatures[ord]
	freq <- freq[ord]
	if(length(freq)==0){
		if(interactive()) warning("All frequency are smaller than min.freq.")
		return(FALSE)
	}
	return(new("motifSig", signatures=signatures, freq=freq, nodelist=nodelist))
}

motifCloud <- function(motifSig, rcpostfix="(RC)", 
layout=c("rectangles", "cloud", "tree"), 
scale=c(6, .5), rot.per=.1,
draw.box=TRUE, draw.freq=TRUE, 
box.col="gray", freq.col="gray",
group.col=NULL, groups=NULL, draw.legend=FALSE)
{
	if (!inherits(motifSig, "motifSig")) 
    stop("motifSig be object of motifSig. You could try\n?motifSignature\nto get a motifSig.")
	signatures=signatures(motifSig)
	freq=frequence(motifSig)
	nodelist=nodelist(motifSig)
	if((!is.null(group.col)) && (!is.null(groups))){
		if(!all(names(table(groups)) %in% names(group.col)))
		stop("some value of groups has no corresponding color in group.col")
		names(groups) <- gsub(rcpostfix,"",names(groups), fixed=T)
		pfmNames <- as.character(unlist(lapply(signatures, function(.ele) unlist(strsplit(.ele@name, ";")))))
		pfmNames <- gsub(rcpostfix,"",pfmNames, fixed=T)
		if(!all(pfmNames %in% names(groups)))
		stop("not all pfms has a given group")
	}
	layout <- match.arg(layout, c("cloud","rectangles", "tree"), several.ok=TRUE)
	layout <- layout[1]
	normedFreq <- freq/max(freq)
	size <- (scale[1]-scale[2])*normedFreq + scale[2]
	grid.pie <- function(x, edges = 200, radius = 0.8, col = NULL, cex=1){
		if (!is.numeric(x) || any(is.na(x) | x < 0)) 
		stop("'x' values must be positive.")
		x <- c(0, cumsum(x)/sum(x))
		dx <- diff(x)
		nx <- length(dx)
		if (is.null(col)) 
		col <- c("white", "lightblue", "mistyrose", "lightcyan", 
				 "lavender", "cornsilk")
		col <- rep(col, length.out = nx)
		twopi <- 2 * pi
		t2xy <- function(t) {
			t2p <- twopi * t
			list(x = radius * cos(t2p)/2 + 0.5, y = radius * sin(t2p)/2 + 0.5)
		}
		for (i in 1L:nx) {
			n <- max(2, floor(edges * dx[i]))
			P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
			grid.polygon(c(P$x, 0.5), c(P$y, 0.5), gp=gpar(col = col[i], fill = col[i]))
		}
	}
	plotSignature <- function(x, y, wid, ht, just, angle, sig, freq, normedFreq){
		pushViewport(viewport(x=x, y=y, width=wid, height=ht, just=just, angle=angle))
		if(draw.box) grid.rect(gp=gpar(col=box.col, lty="dashed", fill="transparent"))
		plotMotifLogoA(sig)
		if(draw.freq) grid.text(label=freq,x=.95,y=.95,gp=gpar(col=freq.col,cex=normedFreq), just=c("right","top"))
		if((!is.null(group.col)) & (!is.null(groups))) {
			sigNames <- unlist(strsplit(sig@name,";"))
			gps <- table(groups[sigNames])
			gp.col <- group.col[names(gps)]
			pushViewport(viewport(x=.1*as.numeric(ht)/as.numeric(wid), y=.9, width=.2*as.numeric(ht)/as.numeric(wid), height=.2))
			grid.pie(gps, col=gp.col, cex=normedFreq)
			popViewport()
		}
		popViewport()
	}
	
	op <- par(mar=c(0,0,0,0))
	on.exit(par(op))
	
	if(layout=="rectangles"){
		plot.new()
		boxesInch <- list()
		minWid <- .Machine$integer.max
		minHt <- .Machine$integer.max
		for(i in 1:length(signatures)){
			ht <- strheight("ACGT", units="inches", cex=size[i])
			wid <- ht*ncol(signatures[[i]]@mat)/2
			boxesInch[[i]] <- new("Rect", width=wid, height=ht)
			if(ht < minHt) minHt <- ht
			if(wid < minWid) minWid <- wid
		}
		pin <- par("pin")
		bin <- new("Rect", width=pin[1], height=pin[2])
		bin.ratio <- bin@height / bin@width
		basic1.wd <- max(boxesInch[[1]]@width, boxesInch[[1]]@height)
		if(bin.ratio < 1) {
			basic1.ht <- basic1.wd * bin.ratio
		}else{
			basic1.ht <- basic1.wd
			basic1.wd <- basic1.ht/bin.ratio
		}
		basic1 <- new("Rect", width=basic1.wd, height=basic1.ht)
		freeRect <- list(basic1)
		bestNodelist <- list()
		maxX <- 0
		maxY <- 0
		for(i in 1:length(signatures)){
			if(length(freeRect)<1){
				if(interactive()) warning("signature could not be fit on page. It will not be plotted.")
				break
			}
			getBestNode <- function(freeRect, boxesInch, i, bin, basic1){
				bestNode <- new("Rect")
				bestX <- basic1@width
				bestY <- basic1@height
				for(j in 1:length(freeRect)){
					if(freeRect[[j]]@width >= boxesInch[[i]]@width && freeRect[[j]]@height >= boxesInch[[i]]@height){
						topSideY <- freeRect[[j]]@y + boxesInch[[i]]@height
						if(topSideY < bestY || (topSideY == bestY && freeRect[[j]]@x < bestX)){
							bestNode@x <- freeRect[[j]]@x
							bestNode@y <- freeRect[[j]]@y
							bestNode@width <- boxesInch[[i]]@width
							bestNode@height <- boxesInch[[i]]@height
							bestY <- topSideY
							bestX <- freeRect[[j]]@x
						}
					}
					if(freeRect[[j]]@width >= boxesInch[[i]]@height && freeRect[[j]]@height >= boxesInch[[i]]@width){
						topSideY <- freeRect[[j]]@y + boxesInch[[i]]@width
						if(topSideY < bestY || (topSideY == bestY && freeRect[[j]]@x < bestX)){
							bestNode@x <- freeRect[[j]]@x
							bestNode@y <- freeRect[[j]]@y
							bestNode@width <- boxesInch[[i]]@height
							bestNode@height <- boxesInch[[i]]@width
							bestY <- topSideY
							bestX <- freeRect[[j]]@x
						}
					}
				}
				if(bestNode@height==0){
					basic2.wd <- min(boxesInch[[i]]@width, boxesInch[[i]]@height)
					basic2.ht <- basic2.wd * bin.ratio
					if(bin.ratio < 1) {
						basic2.ht <- basic2.wd * bin.ratio
					}else{
						basic2.ht <- basic2.wd
						basic2.wd <- basic2.ht/bin.ratio
					}
					basic2 <- new("Rect", width=basic1@width+basic2.wd, height=basic1@height+basic2.ht)
					freeRect <- lapply(freeRect, function(.ele){
									   if(.ele@x+.ele@width==basic1@width) .ele@width <- .ele@width + basic2.wd
									   if(.ele@y+.ele@height==basic1@height) .ele@height <- .ele@height + basic2.ht
									   .ele
									   })
					freeRect <- c(freeRect, new("Rect", x=basic1@width, y=0, width=basic2.wd, height=basic2@height), new("Rect", x=0, y=basic1@height, width=basic2@width, height=basic2.ht))
					basic1 <- basic2
					bestNodeRes <- Recall(freeRect, boxesInch, i, bin, basic1)
					bestNode <- bestNodeRes[["bestNode"]]
					freeRect <- bestNodeRes[["freeRect"]]
					basic1 <- bestNodeRes[["basic1"]]
				}
				list(bestNode=bestNode, freeRect=freeRect, basic1=basic1)
			}
			bestNodeRes <- getBestNode(freeRect, boxesInch, i, bin, basic1)
			bestNode <- bestNodeRes[["bestNode"]]
			freeRect <- bestNodeRes[["freeRect"]]
			basic1 <- bestNodeRes[["basic1"]]
			if(bestNode@height==0){
				if(interactive()) warning("signature could not be fit on page. It will not be plotted.")
				break
			}
			tobeAppend <- list()
			leftNodes <- lapply(freeRect, function(freeNode, usedNode){
								if(usedNode@x >= freeNode@x + freeNode@width || usedNode@x + usedNode@width <= freeNode@x ||
								   usedNode@y >= freeNode@y + freeNode@height || usedNode@y + usedNode@height <= freeNode@y){
								TRUE
								}else{
								if(usedNode@x < freeNode@x + freeNode@width && usedNode@x + usedNode@width > freeNode@x){
## New node at the top side of the used node
								if(usedNode@y > freeNode@y && usedNode@y <freeNode@y + freeNode@height){
								newFreeNode <- freeNode
								newFreeNode@height <- usedNode@y - newFreeNode@y
								tobeAppend <<- c(tobeAppend, newFreeNode)
								}
## New node at the bottom side of the used node
								if(usedNode@y + usedNode@height < freeNode@y + freeNode@height){
								newFreeNode <- freeNode
								newFreeNode@y <- usedNode@y + usedNode@height
								newFreeNode@height <- freeNode@y + freeNode@height - (usedNode@y + usedNode@height)
								tobeAppend <<- c(tobeAppend, newFreeNode)
								}
								}
								if(usedNode@y < freeNode@y + freeNode@height && usedNode@y + usedNode@height > freeNode@y){
## New node at the left side of the used node
								if(usedNode@x > freeNode@x && usedNode@x <freeNode@x + freeNode@width){
								newFreeNode <- freeNode
								newFreeNode@width <- usedNode@x - newFreeNode@x
								tobeAppend <<- c(tobeAppend, newFreeNode)
								}
## New node at the right side of the used node
								if(usedNode@x + usedNode@width < freeNode@x + freeNode@width){
								newFreeNode <- freeNode
								newFreeNode@x <- usedNode@x + usedNode@width
								newFreeNode@width <- freeNode@x + freeNode@width - (usedNode@x + usedNode@width)
								tobeAppend <<- c(tobeAppend, newFreeNode)
								}
								}
								FALSE
								}
								}, bestNode)     
			freeRect <- c(freeRect[unlist(leftNodes)], tobeAppend)
## prune freeRect
			if(length(freeRect)>1){
				leftNodes <- rep(TRUE, length(freeRect))
				for(m in 1:(length(freeRect)-1)){
					if(leftNodes[m]){
						for(n in (m+1):length(freeRect)){
							if(leftNodes[n]){
								if(isContainedIn(freeRect[[m]], freeRect[[n]])){
									leftNodes[m] <- FALSE
									break;
								}
								if(isContainedIn(freeRect[[n]], freeRect[[m]])){
									leftNodes[n] <- FALSE
								}
							}
						}
					}
				}
				freeRect <- freeRect[leftNodes]
			}
			
##push bestNode
			bestNodelist <- c(bestNodelist, bestNode)
			maxX <- max(unlist(lapply(freeRect, function(.ele){.ele@x+.ele@width})))
			maxY <- max(unlist(lapply(freeRect, function(.ele){.ele@y+.ele@height})))
		}
		ratio <- min(bin@width/maxX, bin@height/maxY)
		for(i in 1:length(bestNodelist)){
			bestNode <- bestNodelist[[i]]
## draw logo in bestNode
			if(bestNode@width<bestNode@height){
				plotSignature(x=unit(ratio*bestNode@x, "inches"), 
							  y=unit(ratio*bestNode@y+ratio*bestNode@height, "inches"),
							  wid=unit(ratio*bestNode@height, "inches"), 
							  ht=unit(ratio*bestNode@width, "inches"), 
							  just=c("left","bottom"), 
							  angle=-90,
							  sig=signatures[[i]], 
							  freq=freq[i], 
							  normedFreq=normedFreq[i])
			}else{
				plotSignature(x=unit(ratio*bestNode@x, "inches"), 
							  y=unit(ratio*bestNode@y, "inches"), 
							  wid=unit(ratio*bestNode@width, "inches"), 
							  ht=unit(ratio*bestNode@height, "inches"), 
							  just=c("left","bottom"),
							  angle=0,
							  sig=signatures[[i]], 
							  freq=freq[i], 
							  normedFreq=normedFreq[i])
			}
		}
	}else if(layout=="cloud"){
		plot.new()
		last <- 1
		overlap <- function(x1, y1, sw1, sh1){
			s <- 0
			if(length(boxes)==0) return(FALSE)
			for(i in c(last, 1:length(boxes))){
				bnds <- boxes[[i]]
				x2 <- bnds[1]
				y2 <- bnds[2]
				sw2 <- bnds[3]
				sh2 <- bnds[4]
				if(x1 < x2) overlap <- x1+sw1 > x2-s
				else overlap <- x2 +sw2 > x1-s
				if(y1 <y2) overlap <- overlap && (y1+sh1 > y2-s)
				else overlap <- overlap && (y2 +sh2 > y1-s)
				if(overlap){
					last <- i
					return(TRUE)
				}
			}
			FALSE
		}
		thetaStep <- .1
		rStep <- .05
		plot.window(c(0,1), c(0,1), asp=1)
		boxes <- list()
		
		for(i in 1:length(signatures)){
			sig <- signatures[[i]]
			rotWord <- runif(1)<rot.per
			r <- 0
			theta <- runif(1, 0, 2*pi)
			x1 <- .5
			y1 <- .5
			ht <- strheight("ACGT", cex=size[i])
			wid <- ht * ncol(sig@mat) / 2
			if(rotWord){
				tmp <- ht
				ht <- wid
				wid <- tmp
			}
			isOverlapped <- TRUE
			while(isOverlapped){
				if(!overlap(x1-.5*wid, y1-.5*ht, wid, ht) &&
				   x1-.5*wid>0 && y1-.5*ht>0 &&
				   x1+.5*wid<1 && y1+.5*ht<1){
					if(rotWord){
						plotSignature(x1, y1, ht, wid, "center", 90, sig, freq[i], normedFreq[i])
					}else{
						plotSignature(x1, y1, wid, ht, "center", 0, sig, freq[i], normedFreq[i])
					}
					boxes[[length(boxes)+1]] <- c(x1-.5*wid, y1-.5*ht, wid, ht)
					isOverlapped <- FALSE
				}else{
					if(r>sqrt(.5)){
						if(interactive()) warning("signature could not be fit on page. It will not be plotted.")
						isOverlapped <- FALSE
					}
					theta <- theta+thetaStep
					r <- r+rStep*thetaStep/(2*pi)
					x1 <- .5+r*cos(theta)
					y1 <- .5+r*sin(theta)
				}
			}
		}
	}else{
#trim nodelist
#get all path for each signature
		sigNames <- as.character(unlist(lapply(signatures, function(.ele) .ele@name)))
		sigNames <- gsub(rcpostfix, "", sigNames, fixed=T)
		paths <- lapply(sigNames, function(.ele, nodelist){
						pt <- c()
						getPath <- function(nodelist, nodename, pt){
						parent <- getParentNode(nodelist, nodename)
						if(!is.null(parent[1])){
						pt <- c(pt, parent[1])
						pt <- Recall(nodelist, parent[1], pt)
						}
						pt
						}
						pt <- getPath(nodelist, .ele, pt)
						pt
						}, nodelist)
#leave the nodes in the paths
		paths <- unique(unlist(paths))
		nodelist <- nodelist[paths]
		
		last <- 1
		overlap <- function(x1, y1, sw1, sh1){
			s <- 0
			if(length(boxes)==0) return(FALSE)
			for(i in c(last, 1:length(boxes))){
				bnds <- boxes[[i]]
				x2 <- bnds[1]
				y2 <- bnds[2]
				sw2 <- bnds[3]
				sh2 <- bnds[4]
				if(x1 < x2) overlap <- x1+sw1 > x2-s
				else overlap <- x2 +sw2 > x1-s
				if(y1 <y2) overlap <- overlap && (y1+sh1 > y2-s)
				else overlap <- overlap && (y2 +sh2 > y1-s)
				if(overlap){
					last <- i
					return(TRUE)
				}
			}
			FALSE
		}
		boxes <- list()
		step <- .1
		plot.default(0, 0, type = "n", asp = 1, xlab = "", ylab = "", 
					 xaxt = "n", yaxt = "n", xlim = c(0, 1), ylim = c(0, 1), 
					 xaxs = "i", yaxs = "i", frame.plot = FALSE)
		getRdist <- function(nodelist, Rdist=c(), nodename="Root", distp=0){
			currNode <- nodelist[[nodename]]
			distl <- distp + currNode@distl
			distr <- distp + currNode@distr
			currDist <- c(distl, distr)
			names(currDist) <- c(currNode@left, currNode@right)
			Rdist <- c(Rdist, currDist)
			if(!is.null(nodelist[[currNode@left]])) Rdist <- getRdist(nodelist, Rdist, currNode@left, distp=distl)
			if(!is.null(nodelist[[currNode@right]])) Rdist <- getRdist(nodelist, Rdist, currNode@right, distp=distr)
			Rdist
		}
		Rdist <- getRdist(nodelist)
		longestRdist <- 4*max(Rdist)
		plotUnrootedTree <- function(signatures, nodelist, positions, nodename="Root", AXIS=0, ANGLE=2*pi, orignalX=0.5, orignalY=0.5, longestRdist=0){
			currNode <- nodelist[[nodename]]
			RdistTotal <- 0
			if(currNode@left %in% c(paths, sigNames))
			RdistTotal <- currNode@distl
			if(currNode@right %in% c(paths, sigNames)) 
			RdistTotal <- RdistTotal + currNode@distr
#plot left node
			start <- AXIS - ANGLE/2
			alpha <- 0
			if(currNode@left %in% c(paths, sigNames)){
				alpha <- ANGLE*currNode@distl/RdistTotal
				beta <- start + alpha/2
				x <- currNode@distl/longestRdist * cos(beta) + orignalX
				y <- currNode@distl/longestRdist * sin(beta) + orignalY
				if(!is.null(nodelist[[currNode@left]])) {
					segments(orignalX, orignalY, x, y)
					positions <- plotUnrootedTree(signatures, nodelist, positions, nodename=currNode@left, AXIS=beta, ANGLE=alpha, orignalX=x, orignalY=y, longestRdist=longestRdist)
				}else{
					pfmlid <- getPFMid(signatures, currNode@left)
					if(length(pfmlid)>0){
						ht <- strheight("ACGT", cex=size[pfmlid])
						wid <- ht * ncol(signatures[[pfmlid]]@mat) / 2
						segments(orignalX, orignalY, x, y)
						cp <- new("Rect", x=x, y=y, width=wid, height=ht)
						positions <- c(positions, new("pos", box=cp, beta=beta, sig=signatures[[pfmlid]], freq=freq[[pfmlid]], norm=normedFreq[pfmlid]))
					}
				}
			}
#plot right node
			if(currNode@right %in% c(paths, sigNames)){
				start <- start + alpha
				alpha <- ANGLE*currNode@distr/RdistTotal
				beta <- start + alpha/2
				x <- currNode@distr/longestRdist * cos(beta) + orignalX
				y <- currNode@distr/longestRdist * sin(beta) + orignalY
				if(!is.null(nodelist[[currNode@right]])) {
					segments(orignalX, orignalY, x, y)
					positions <- plotUnrootedTree(signatures, nodelist, positions, nodename=currNode@right, AXIS=beta, ANGLE=alpha, orignalX=x, orignalY=y, longestRdist=longestRdist)
				}else{
					pfmrid <- getPFMid(signatures, currNode@right)
					if(length(pfmrid)>0){
						ht <- strheight("ACGT", cex=size[pfmrid])
						wid <- ht * ncol(signatures[[pfmrid]]@mat) / 2
						segments(orignalX, orignalY, x, y)
						cp <- new("Rect", x=x, y=y, width=wid, height=ht)
						positions <- c(positions, new("pos", box=cp, beta=beta, sig=signatures[[pfmrid]], freq=freq[[pfmrid]], norm=normedFreq[pfmrid]))
					}
				}
			}
			positions
		}
		positions <- plotUnrootedTree(signatures, nodelist, positions=list(), longestRdist=longestRdist)
		positions <- positions[order(unlist(lapply(positions, function(.ele) .ele@freq)), decreasing=T)]
		for(i in 1:length(positions)){
			isOverlapped <- TRUE
			thisMotif <- positions[[i]]
			x <- thisMotif@box@x
			y <- thisMotif@box@y
			wid <- thisMotif@box@width
			ht <- thisMotif@box@height
			beta <- thisMotif@beta
			x1 <- x
			y1 <- y
			while(isOverlapped){
				if(!overlap(x1-.5*wid, y1-.5*ht, wid, ht)){
					segments(x,y,x1,y1,lty="dashed", col="gray")
#plot signature
					plotSignature(x1, y1, wid, ht, "center", 0, thisMotif@sig, thisMotif@freq, thisMotif@norm)
					boxes[[length(boxes)+1]] <- c(x1-.5*wid, y1-.5*ht, wid, ht)
					isOverlapped <- FALSE
				}else{
					if((x1>=1-.5*wid || x1<=.5*wid) && (y1<=.5*ht || y1>=1-.5*ht)){
						if(interactive()) warning("signature could not be fit on page. It will not be plotted.")
						isOverlapped <- FALSE
					}
					r <- sqrt(wid*wid + ht*ht)*step
					if(((beta %% pi) - pi/2) < 0.0001) beta <- beta + pi/90
					x1 <- x1+r*cos(beta)
					y1 <- y1+r*sin(beta)
					if(x1-.5*wid<0) x1 <- .5*wid
					if(x1+.5*wid>1) x1 <-1-.5*wid
					if(y1-.5*ht<0) y1 <- .5*ht
					if(y1+.5*ht>1) y1 <-1-.5*ht
				}
			}
		}
	}
#plot group legend
	if((!is.null(group.col)) && (!is.null(groups)) && draw.legend){
		legend("topright", names(group.col), fill=group.col, cex=0.5)
	}
	return(invisible())
}