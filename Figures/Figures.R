require(Cairo)

ppi<-300
png("Fig1_GEBCO_Balt.png", width = 3.5*ppi, height = 3*ppi, res=ppi,type ="cairo",bg = "transparent")
par(mar=c(2,2,2,2))
plot(bathy, xlim=c(13, 22), ylim=c(54, 59),
     col = brewer.pal(9,"Blues"), main="Bathymetry")
map("worldHires", fill=TRUE, col="gray90", border="gray70", add=TRUE)
plot(hauls, add=TRUE)
mtext(text = "Depth (m)", side = 4, line = 0.2, cex=0.7)
box()
dev.off()

onetwo <- function(dat){
  return(paste(dat[1], dat[2])) 
}
shortnam <- unlist(lapply(strsplit(levels(habitats$Grouped), " "),onetwo))
png("Fig2_Seabed.png", width = 3*ppi, height = 3*ppi, res=ppi,type ="cairo",bg = "transparent")
par(mar=c(2,2,2,1))
plot(habitats, col=brewer.pal(9, "Set3")[habitats$Grouped], 
     border=FALSE, main="Seabed habitat", xlim=c(13, 22), ylim=c(54, 59))
plot(hauls, add=TRUE, cex=0.5)
map("worldHires", col="gray90", add=TRUE)
legend("topleft", legend = shortnam, bty="n",
       fill = brewer.pal(9, "Set3"), cex=0.6)
box()
axis(1)
axis(2)
dev.off()

png("Fig3_GlobColour_NS.png", width = 3.9*ppi, height = 3*ppi, res=ppi,type ="cairo",bg = "transparent")
par(mar=c(2,2,3,2))
plot(GColor072015, xlim=c(4.5,9.5), ylim=c(53.2,55.7), col=brewer.pal(9, "Greens"),
     main="GlobColour", sub="CHl concentration (mg/m3)")
map("worldHires", fill=TRUE, col="gray90", border="gray70", add=TRUE)
plot(ICESrect, add=TRUE)
text(ICESrect, label=ICESrect$ICESNAME, cex=0.8)
box()
mtext(text = "CHl concentration (mg/m3)", side = 4, line = 0.2, cex=0.7)
dev.off()

png("Fig4_SufTemp_NS.png", width = 3.9*ppi, height = 3*ppi, res=ppi,type ="cairo",bg = "transparent")
plot(temp072015, y=1, main="Surface temperature in July 2015", 
     col= brewer.pal(9,"YlOrRd"))
plot(hauls, add=TRUE, cex=0.5)
map("worldHires", fill=TRUE, col="gray90", border="gray70",
    add=TRUE)
dev.off()

png("Fig5_DepthTemp.png", width = 3.5*ppi, height = 3.5*ppi, res=ppi,bg = "transparent")
par(mar=c(4,4,2,0.5))
pal <- brewer.pal(9, "YlOrRd")
coldp10 <- colscale(depth10, pal)
plot(haul.temp[1,], -depth , type="l", xlim=range(haul.temp, na.rm=TRUE),
     xlab="Temperature", ylab="Depth", las=1, col=coldp10$col[1], 
     main="Thermocline")
for (i in 2:nrow(haul.temp)){
  lines(haul.temp[i,],-depth, col=coldp10$col[i])
}
abline(v=10, lty=3)
dev.off()

png("Fig6_DepthTemp.png", width = 3*ppi, height = 3*ppi, res=ppi,bg = "transparent")
#par(mar=c(2,2,2,1))
#size <- 2*log(hauls$totCPUE)/log(max(hauls$totCPUE))
map("worldHires",  xlim=c(14,22), ylim=c(54,58.5), 
    fill=TRUE, col="gray90", border="gray70", mar = c(2,2,2,1))
plot(hauls, col=coldp10$col, pch=16, add=TRUE)
box()
axis(1)
axis(2)
title("Depth at 10°C")
add.colscale(coldp10$br, pal,posi="topleft", lab="Depth (m)", rd=0, ratio = 0.3, cex=0.6)
dev.off()

png("Fig6_SST.png", width = 3.5*ppi, height = 3.8*ppi, res=ppi,bg = "transparent")
par(mar=c(2,2,2,0.5))
plot(temp072015, y=1, main="Surface temperature", 
     col= brewer.pal(9,"YlOrRd"))
plot(hauls, add=TRUE, cex=0.4)
map("worldHires", fill=TRUE, col="gray90", border="gray70",
    add=TRUE)
dev.off()

pal <- brewer.pal(9, "Greens")
colNPP <- colscale(maxNPP, pal)
png("Fig7_ChlSea.png", width = 4*ppi, height = 3*ppi, res=ppi,bg = "transparent")
par(mar=c(2,4,2,0.5))
plot(meanchl[1,], ylim=range(meanchl), lwd=2, xaxt="n", main="Seasonal Primary Production",
     type="l", col=colNPP$col[1], ylab="CHL concentration (mg/m3)")
axis(side = 1, at = 1:11, labels = month.abb[1:11])
for (i in 2:nrow(meanchl)){
  lines(meanchl[i,], ylim=range(meanchl), lwd=2, 
        type="l", col=colNPP$col[i])
}
dev.off()


maxNPP<- apply(meanchl, 1, max)

maxSea<- apply(meanchl, 1, which.max)
maxSea.abb <- month.abb[maxSea]

pal <- brewer.pal(9, "Greens")
colNPP <- colscale(maxNPP, pal)
map("worldHires",  xlim=c(4,10), ylim=c(53,56), 
    fill=TRUE, col="gray90", border="gray70")
plot(ICESrect, col=colNPP$col, add=TRUE)
text(ICESrect, maxSea.abb)
box()
axis(1)
axis(2)
title("Seasonal pic of primary production")
add.colscale(colNPP$br, pal,posi="bottomright", lab="chl (mg/m3)", rd=0)