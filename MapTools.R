#Home-made functions to create color scale and plot them
#R. Frelat - 29.06.2017

#colscale : compute the scale of colors uniformly from the range of values
# dat : vector of data
# pal : color palette
colscale <- function(dat, pal, na.rm=TRUE){
  br <- seq(min(dat, na.rm = na.rm), max(dat, na.rm = na.rm), length.out = length(pal)+1)
  col <- pal[cut(dat, include.lowest = TRUE, breaks=br)]
  res <- list("br"=br, "col"=col)
  return(res)
}

#add.colscale : add the scale of colors in a plot
# dat : vector of data
# pal : color palette
add.colscale <- function(br, pal,posi = "topleft", rd=1, las=1, lab="", cex=0.8, ratio = 0.2, inset = 0.01){
  add.scatter(plot.colscale(br, pal, rd=rd, las=las, lab=lab, cex=cex, posi=posi), posi=posi, ratio = ratio, inset = inset)
}

#plot.scale : function to plot the scale in a separate plot
plot.scale <- function(br, pal, rd=1, n=4, lab="", cex=0.8, ...){
  plot(0, xlim=c(0,1), ylim=c(0.1,length(pal)-0.1), type="n", 
       xaxt="n", yaxt="n", bty="n", xlab="", ylab="", ...)
  for (i in 1:length(pal)){
    rect(0,i-1,1,i,col = pal[i])
  }
  showAx <- round(seq(0,length(pal), length.out=n))
  axis(2, at=showAx, labels = round(br,rd)[showAx+1])
  mtext(lab, 3, line = 0.5, adj = 1, cex=cex)
}

#plot.colscale : function to be coupled with add.scatter
plot.colscale <- function(br, pal, rd=1, las=1, posi="left", lab="", cex=0.8){
  opar=par("mar","xaxt","yaxt","plt")
  on.exit(par(opar))
  axp <- ifelse(posi == "bottomright" || posi == "topright", 2, 4)
  adjp <- ifelse(posi == "bottomright" || posi == "topright", 1, 0)
  par(mar=c(0,2,2,2),plt=par("plt"), las=las)
  plot(0, xlim=c(0,1), ylim=c(0.1,length(pal)-0.1), type="n", 
       bty="n", xlab="", ylab="", xaxt="n",yaxt="n")
  for (i in 1:length(pal)){
    rect(0,i-1,1,i,col = pal[i])
  }
  showAx <- round(seq(0,length(pal), length.out=4))
  axis(axp, at=showAx, labels = round(br,rd)[showAx+1], cex.axis=cex)
  mtext(lab, 3, line = 0.5, adj = adjp, cex=cex)
}

#add.scatter : function modified from ade4 package
add.scatter <- function (func, posi = c("bottomleft", "bottomright", "topleft", 
                         "topright"), ratio = 0.2, inset = 0.01, bg.col="white") {
  if (tolower(posi[1]) == "none") 
    return()
  if (ratio > 0.99) 
    ratio <- 0.99
  if (ratio < 0) 
    ratio <- 0.2
  if (length(inset) == 2) {
    inset.x <- inset[1]
    inset.y <- inset[2]
  }
  else {
    inset.x <- inset[1]
    inset.y <- inset[1]
  }
  inset[inset < 0] <- 0
  plotreg0 <- par("plt")
  plotreg <- plotreg0 + c(inset.x, -inset.x, inset.y, -inset.y)
  on.exit(par(plt = plotreg0))
  posi <- tolower(posi[1])
  if (posi == "bottomleft" || posi == "bottom") {
    x1 <- plotreg[1]
    y1 <- plotreg[3]
  }
  else if (posi == "topleft" || posi == "top") {
    x1 <- plotreg[1]
    y1 <- plotreg[4] - ratio
  }
  else if (posi == "bottomright") {
    x1 <- plotreg[2] - ratio
    y1 <- plotreg[3]
  }
  else if (posi == "topright") {
    x1 <- plotreg[2] - ratio
    y1 <- plotreg[4] - ratio
  }
  else stop("Unknown position required")
  x2 <- x1 + ratio
  y2 <- y1 + ratio
  # if (posi == "bottomright" || posi == "topright") {
  #   par(plt = c(x1+(ratio*0.5), x2, y1, y2), new = TRUE)
  # } else {
  #   par(plt = c(x1, x1 + ratio/2, y1, y2), new = TRUE)
  # }
  # plot.new()
  # polygon(c(-0.1, 1.1, 1.1, -0.1), c(-0.1, -0.1, 1.1, 1.1), 
  #         border = NA, col = bg.col)
  if (posi == "bottomright" || posi == "topright") {
    par(plt = c(x1+(ratio*0.8), x2, y1+ (ratio/10), y2- (ratio/4)), new = TRUE)
  } else {
    par(plt = c(x1, x1 + ratio*0.2, y1+ (ratio/10), y2- (ratio/4)), new = TRUE)
    # par(plt = c(x1+(ratio*0.4), x1 + ratio*0.6, y1+ (ratio/10), y2- (ratio/4)), new = TRUE)
  }
  eval(func)
  return(invisible(match.call()))
}