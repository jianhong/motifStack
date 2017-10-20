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