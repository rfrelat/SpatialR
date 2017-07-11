## Workshop : Spatial analysis with R
# Hands-on multidimensional examples
# By Romain Frelat, 10 July 2017

## A. Introduction: ----------------------------
# To load and process GIS data
require(sp)
require(rgdal)
require(raster)
require(ncdf4)
#To make nicer looking maps
require(maps) 
require(mapdata)
require(RColorBrewer)
source("MapTools.R")

## B. Depth profile of the Baltic Sea -----------

# B.1 Load and vizualize the data
dir <- "Data/CMEMS_SMHI_PHYS_reanalysis_201507.nc"
temp072015 <- brick(dir, varname="temp", lvar=4)

#Access information from the brick object
proj4string(temp072015)
dim(temp072015)
res(temp072015)

#Depth layers
depth <- temp072015@z[[1]]
#See the depth layers
depth

# Load the hauls shapefile
dir <- "Data/Hauls_BITSTVLQ1_2016.shp"
name <- "Hauls_BITSTVLQ1_2016"
hauls <- readOGR(dir,name)

# Visualize
plot(temp072015, y=1, main="Surface temperature in July 2015", 
     col= brewer.pal(9,"YlOrRd"))
plot(hauls, add=TRUE, cex=0.5)
map("worldHires", col="grey90", border="grey50", 
    fill=TRUE, add=TRUE)

#Your turn : ....................................

# B.2 Extract values
haul_temp <- extract(temp072015, hauls)
haul_temp <- extract(temp072015, hauls, fun=mean)
dim(haul_temp) #158 hauls x 50 depth layer
colnames(haul_temp) <- depth

#Remove depth layer deeper than 120m
haul_temp <- haul_temp[,depth<=120]
depth <- depth[depth<=120]

# B.3 Get surface and bottom temperature
#Surface temperature : first column
sst <- haul_temp[,1] 

#Identify the bottom for each haul
bottom <- apply(!is.na(haul_temp), 1, sum)
#Bottom temperature
sbt <- haul_temp[cbind(1:nrow(haul_temp),bottom)]

#Choose the color palette
pal <- brewer.pal(9, "YlOrRd")
# Assign a color to each haul according to sst
col_sst <- colscale(sst, pal)
#Create a map
#Choose the color palette
pal <- brewer.pal(9, "YlOrRd")
# Assign a color to each haul according to sst
col_sst <- colscale(sst, pal)
#Create a map
plot(hauls, col=col_sst$col, pch=16, axes=TRUE,
     main="Surface temperature in July 2015")
map("worldHires", col="grey90", border="grey50", 
    fill=TRUE, add=TRUE)
#Add the color scale
add.colscale(col_sst$br, pal, posi="topleft", 
             lab="Temperature", cex=0.7)

#Your turn : ....................................

# B.4 Visualize the thermocline
plot(haul_temp[1,], -depth , type="l", , col="grey30",
     xlim=range(haul_temp, na.rm=TRUE), xlab="Temperature",
     ylab="Depth", las=1)
for (i in 2:nrow(haul_temp)){
  lines(haul_temp[i,],-depth, col="grey30")
}

#Detect the depth when the temperature pass the 10 degrees threshold
depth10 <- depth[apply(haul_temp>10, 1, sum, na.rm=TRUE)]
#remove the hauls which never reach lower temperature than 10 degrees
depth10[sbt>10] <- NA

#Visualize the distribution of the depth at 10 degrees
boxplot(depth10, ylab="Depth (m)", main="Depth at 10 degrees")

#Your turn : ....................................

## C Seasonality of primary production in the North Sea

# C1. Load and vizualize the data
#list all the files in GlobColour folder finishing by "00.nc"
file.names<-list.files("Data/GlobColour", pattern="00.nc$", 
                       full.names = TRUE)
length(file.names)
# Getting the date of each file from the file name
# The date is between position 21 and 26.
time<-substr(file.names, 21,26)
time

# Load a stack object
GColor2015<-stack(file.names, varname="CHL1_mean")
names(GColor2015) <- time

#Access information from the stack object
proj4string(GColor2015)
dim(GColor2015)
res(GColor2015)

# Load the ICES rectangle
dir <- "Data/ICESrect_GermanBight.shp"
name <- "ICESrect_GermanBight"
ICESrect <- readOGR(dir,name)

# Visualize
#Define manually the color breaks
brk <- c(0,1,3,6,10,15,20,25,30,35)
#Define the color palette
pal <- brewer.pal(9, "BuGn")
#Plot the 11 months
par(mfrow=c(3,4), mar=c(2,2,2,1))
for (i in 1:length(time)){
  #show the primary production of time i
  image(GColor2015, y=i, main="", xlim=c(4,10),  
        ylim=c(53, 56), col=pal, breaks=brk)
  #add the country borders
  map("worldHires", col="grey90", border="grey50", 
      fill=TRUE, add=TRUE)
  #add the ICES rectangle
  plot(ICESrect, add=TRUE)
  #add the title - the time
  title(time[i], adj=1)
}
#Add the color scale in a separate plot
par(mar=c(2,8,4,4))
plot.scale(brk, pal = pal, lab="Chl concentration\n(mg/m3)")

# C2. Extract monthly primary production in ICES rectangle

#compute the sum per polygon
ICESGColor_sum <- extract(GColor2015, ICESrect, 
                          fun=sum, na.rm=TRUE)
#get the number of pixel with non-na values
lena <- function(dat, ...){return(sum(!is.na(dat)))}
ICESGColor_length <- extract(GColor2015, ICESrect, fun=lena)
#calculate the correct average
ICESGColor_mean <- ICESGColor_sum/apply(ICESGColor_length,1,max)
dim(ICESGColor_mean)

# C3. Visualize the seasonal pattern
#Layout to match the ICES spatial distribution
layout(matrix(c(4:1, 8:5, 12:9, 16:13), byrow = FALSE, ncol=4))
#Plot the seasonal primary production
par(mar=c(2,2,2,0))
for (i in 1:nrow(ICESGColor_mean)){
  plot(ICESGColor_mean[i,], ylim=range(ICESGColor_mean), lwd=2, 
       type="l", main=ICESrect$ICESNAME[i])
}

#Your turn : ....................................


## End !!