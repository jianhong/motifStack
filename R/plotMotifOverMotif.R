## compare two motifs
## require two motifs must be exact same width and aligned
plotMotifOverMotif <- function(motif, backgroundMotif, bgNoise=NA, 
                               font="Helvetica-Bold", textgp=gpar()){
    if(!inherits(motif, c("pcm", "pfm")))
        stop("motif must be an object of pcm or pfm")
    if(!inherits(backgroundMotif, c("pcm", "pfm")))
        stop("backgroundMotif must be an object of pcm or pfm")
    if(class(motif)!="pfm") motif <- pcm2pfm(motif)
    if(class(backgroundMotif)!="pfm") backgroundMotif <- pcm2pfm(backgroundMotif)
    
    w1 <- ncol(motif$mat)
    w2 <- ncol(backgroundMotif$mat)
    if(w1!=w2) stop("motif and backgroundMotif must be exact same width and aligned")
    if(!identical(motif$alphabet, backgroundMotif$alphabet))
        stop("the type of motif and backgroundMotif is different.")
    if(!identical(motif$background, backgroundMotif$background))
        stop("the backgrounds of motif and backgroundMotif are not identical.")
    
    pfm <- motif$mat
    background <- backgroundMotif$mat
    rdirichlet <- function (n, alpha){
        l <- length(alpha)
        x <- matrix(rgamma(l * n, alpha), ncol = n, byrow = TRUE)
        sm <- x %*% rep(1, n)
        t(x/as.vector(sm))
    }
    if(!is.na(bgNoise)){
        if(bgNoise<0 | bgNoise>1) stop("bgNoise must be number in (0, 1)")
        background <- (1-bgNoise)*background + 
            bgNoise*rdirichlet(nrow(background), rep(1, ncol(background)))
    }
    #plotMotifLogo(pfm, motifName, background, ...)
    #plot
    dat <- pfm - background
    gscmd <- Sys.getenv("R_GSCMD")
    npos <- ncol(dat)
    ncha <- nrow(dat)
    colset <- motif$color
    rname <- rownames(dat)
    key<-paste("x", ncha, font, paste(colset[rname], collapse=""), 
               paste(rname, collapse=""), sep="_")
    symbolsCache <- 
        if(exists("tmp_motifStack_symbolsCache", where=".GlobalEnv")){
            get("tmp_motifStack_symbolsCache", pos=".GlobalEnv")
        } else {
            list()
        }
    if(!is.null(symbolsCache[[key]])){
        symbols<-symbolsCache[[key]]
    } else {
        symbols<-motifStack:::coloredSymbols(ncha, font, colset[rname], rname)
        symbolsCache[[key]]<-symbols
        assign("tmp_motifStack_symbolsCache", symbolsCache, pos=".GlobalEnv")
    }
    
    ##check font height to fix the number size
    plot.new()
    dw <- 1/(npos+2)
    line1 <- as.numeric(convertUnit(unit(1, "strwidth", "W"), "npc"))
    pin <- dev.size("in") #pin <- par("pin")
    cex <- dw/line1
    if(length(textgp)==0){
        cex1 <- ifelse(pin[1L] > pin[2L],
                       cex * pin[2L]/pin[1L],
                       cex)
        textgp <- gpar(cex=.8 * cex1)
    }
    tickgp <- textgp
    tickgp$cex <- cex/3
    x0 <- dw ## max 4 chars
    x1 <- 4/5 * dw
    x2 <- 6/5 * dw
    ##set ylim
    datN <- apply(dat, 2, function(.col) sum(.col[.col<0]))
    datP <- apply(dat, 2, function(.col) sum(.col[.col>0]))
    ylim <- c((as.integer(min(datN)/0.05)-1)*0.05, (as.integer(max(datP)/0.05)+1)*0.05)
    remap <- function(x){
        (ylim[1] - x)/(ylim[1] - ylim[2])/(1+dw)
    }
    reheight <- function(x){
        abs(x/(ylim[1] - ylim[2])/(1+dw))
    }
    
    ##draw axis
    x.pos <- 0
    grid.text(0, x0, remap(0)+dw/2, just=c(.5, .5), gp=textgp)
    grid.lines(x=x0, y=c(remap(0)+dw, 1), arrow=arrow(length=unit(0.02, "npc")), gp=gpar(lwd=tickgp$cex))
    grid.lines(x=x0, y=c(remap(0), 0), arrow=arrow(length=unit(0.02, "npc")), gp=gpar(lwd=tickgp$cex))
    ##draw tick
    tick <- 0.1
    times <- 100
    for(i in c(as.integer(min(datN)/tick):(-1), 1:as.integer(max(datP)/tick))){
        grid.lines(x=c(x2, x0), y=remap(i*tick)+dw/2, gp=gpar(lwd=tickgp$cex))
        grid.text(times*tick*i, x=x1, remap(i*tick)+dw/2, just=c(1, .5), gp=tickgp)
    }
    
    x.pos <- x0 + dw/2
    for(j in 1:npos){
        heights <- dat[, j]
        id <- order(heights)
        id1 <- order(heights)
        y.pos <- remap(sum(heights[heights<0]))
        flag <- 0
        if(length(heights)>0){
            for(i in 1:length(heights)){
                h <- reheight(heights[id1[i]])
                if(heights[id1[i]]>0) flag <- flag+1
                if(flag==1) {
                        grid.text(j, x.pos+dw/2, y.pos+dw/2, just=c(.5, .5), gp=textgp)
                    y.pos <- y.pos + dw
                    flag <- flag+1
                }
                if(h>0) {
                    grid.draw(grImport::pictureGrob(symbols[[id[i]]],
                                                    x.pos,y.pos,dw,h,
                                                    just=c(0,0),distort=TRUE))
                    y.pos<-y.pos+h
                }
            }
            
        }
        if(flag==0) grid.text(j, x.pos+dw/2, y.pos+dw/2, just=c(.5, .5), gp=textgp)
        x.pos<-x.pos+dw
    }
    return(invisible(dat))
}