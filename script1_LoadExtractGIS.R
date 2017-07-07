## Workshop : Spatial analysis with R
# How to load, extract and analyse spatial data in R ?
# By Romain Frelat, 10 July 2017

## A. Getting ready: ----------------------------
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

## B. Load spatial data -------------------------

# B.1 Load a raster

# Bathymetry of the baltic Sea
#Set the directory of the file
dir<-"Data/GEBCO2014_Subset_30Sec.tif" 
#Load the raster
bathy <- raster(dir) 

# Information about the loaded raster
#Projection system
proj4string(bathy)
#Extent of the dataset
bbox(bathy)
#Dimension and resolution
dim(bathy)
res(bathy)

# Visualize 
#Visualize the raster
plot(bathy, maxpixels = 20000)
#Add the country borders
map("worldHires", col="grey90", border="grey50", 
    fill=TRUE, add=TRUE)

#Your turn : ....................................

# B.2 Load a vector 

# Bottom trawl survey in the baltic Sea
dir <- "Data/Hauls_BITSTVLQ1_2016.shp"
name <- "Hauls_BITSTVLQ1_2016"
hauls <- readOGR(dir,name)

# Projection and extent
#See the projection
proj4string(hauls)
#See the extent
bbox(hauls)
#Dimension of the attribute table
dim(hauls)
#Variables in the attribute table
names(hauls@data)

#Identify which species is the most abundant:
species <- c("Herring", "Cod", "Flounder", "Plaice", "Sprat")
abu <- hauls@data[,species]
dom_sp <- as.factor(species[apply(abu, 1, which.max)])
table(dom_sp)

# Visualize the shapefile
plot(hauls)
#Add a background of countries border
map("worldHires", col="grey90", border="grey50", fill=TRUE, add=TRUE)
#Add a box and axis
box()
axis(1)
axis(2)

# log transform CPUE
logCPUE <- log(hauls@data$totCPUE)

#Size of the dot, proportional to the log of total catch
size <- 2*logCPUE/max(logCPUE) #ratio between 0 and 2 of the total catch

#Set the palette of colors for the 4 dominant species
pal <- brewer.pal(4, "Set2")

#Plot the hauls with 'cex' telling the size of the dot, and 'col' the color
plot(hauls, pch=16, cex=size, col=pal[dom_sp])

#Add a background, a legend and borders
map("worldHires", col="grey90", border="grey50", fill=TRUE, add=TRUE)
legend("topleft", legend = levels(dom_sp), col = pal,
       pch = 16, title="dom. species")
box()
axis(1)
axis(2)

#Your turn : ....................................

## C. Crossing points and raster ----------------

# Visually overlay raster and points
plot(bathy, xlim=c(13, 22), ylim=c(54, 59),
     col = brewer.pal(9,"Blues"), main="Bathymetry and hauls")
plot(hauls, add=TRUE)
map("worldHires", col="grey90", border="grey50", 
    fill=TRUE, add=TRUE)

# Extract information from a raster
hauls_depth <- extract(bathy, hauls)
length(hauls_depth)

# Visualize extracted information
plot(hauls$Depth, hauls_depth, xlab="Depth from BITS",
     ylab="Depth from GEBCO")
abline(a=0, b=1) #add the y=x line

#Your turn : ....................................

# D. Crossing points with polygons --------------

# Seabed habitats
dir <- "Data/Baltic_Habitats_EUSEaMap2011.shp"
name <- "Baltic_Habitats_EUSEaMap2011"
habitats <- readOGR(dir,name)

#For Mac users only
require(Cairo)
X11(type="cairo")

#Your turn : ....................................

# Crossing points with polygons
hauls_hab <- over(hauls, habitats) 
dim(hauls_hab)


## E. Cross polygons with raster

# Load GlobColour raster
dir <- "Data/GlobColour/L3m_20150701-20150731__728672765_1_GSM-MODVIR_CHL1_MO_00.nc"
GColor072015 <- raster(dir, varname="CHL1_mean")

# Load ICES rectangle
dir <- "Data/ICESrect_GermanBight.shp"
name <- "ICESrect_GermanBight"
ICESrect <- readOGR(dir,name)

#Your turn : ....................................

# Crossing polygons with raster

# Extract mean values per polygon
ICESGColor_mean <- extract(GColor072015, ICESrect, fun=mean, na.rm=TRUE)
dim(ICESGColor_mean)

# Visualize the extracted information
#Choose the palette colors
pal <- brewer.pal(9, "Greens")
#Create the color scale
col_mean <- colscale(ICESGColor_mean, pal)
#Map the rectangle according to the color scale
plot(ICESrect, col=col_mean$col)
#Add country border and axis
map("worldHires", col="black", add=TRUE)
box()
axis(1)
axis(2)
#Add a color scale on the map
add.colscale(col_mean$br, pal,posi="topleft", lab="Chl (mg/m3)")

#Extract all the pixel values per polygon
ICESGColor <- extract(GColor072015, ICESrect)
typeof(ICESGColor)
length(ICESGColor)
names(ICESGColor) <- ICESrect$ICESNAME

# Visualize the variation within each polygon
boxplot(ICESGColor, las=2)
#add a line with the previously extracted mean
lines(ICESGColor_mean, lwd=2, col="red")

## End !!