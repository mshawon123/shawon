# Spatial Analysis of landset ARD Data
# Project-1
# Mostafizur Shawon
# Decade 2014-2023

# Step:1 load libraries
library(terra)
library(gdata)
library(corrplot)
library(sf)
library(pracma)
library(dplyr)
library(ggplot2)
library(cluster)
library(stats)
library(factoextra)


#Step:2 Set working directory
path <- "D:\\Remote sensing\\Project_1\\Dataset_2023\\Bands"
setwd(path)

#Step 3: read all the bands
bands <- c('aerosolx','bluex','redx','greenx','NIRx','SWR1x','SWR2x')
Nbands<- length(bands)

#Step 4: Get ready to read the hardin vector file
fname <- 'Hardin.gpkg' 
hardin <- vect(fname)
plot (hardin)

#Step 5: Read the raster and crop 
for (i in seq(1,Nbands,1)){
  fname <- paste('Band',i,'.tif',sep="")
  x <- rast(fname)
  x <- crop(x,hardin)
  plot(x)
  mv('x',bands[i])
}
#Step 6: Create a rasterstack with all data
hardinlsat <- c(aerosolx,bluex,greenx,redx,NIRx,SWR1x,SWR2x)
writeRaster (hardinlsat, "hardinallbands.TIF" ,
             filetype = "GTiff", overwrite = TRUE)

# Step 7: Create a Boxplot
bandnames <- c( 'Aerosol','Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2')
f <- boxplot(hardinlsat,axes=FALSE,outline=F,ylab='Value',notch=F)
axis(1, at=1: (Nbands),labels=bandnames, las=2)
axis (2)
title('Digital numbers for Various Bands in Hardin')
grid()
box()

###################################################################################
#Step 8: Write the band as data frame for additional analysis
banddf <- as.data.frame(hardinlsat, xy = T)
head(banddf)
banddf <- na.omit(banddf)
colnames(banddf) <- c('X', 'Y', bandnames)
head(banddf)


#Step 9 :correlation between bands
corbands <- cor(banddf[,3:8], method = 'pearson') #spearman
corrplot(corbands, method = 'number', type = 'lower', diag = F)


#step 10: Clamp values for contrast enhancement #histogram equalization
LL <- 0.05
UL <- 0.95
bandnamec <- c()
for (i in seq(1, Nbands, 1)){
  bname <- paste(bands[i], 'c', sep = "")
  xmin <- quantile(banddf[, (i+2)], LL, na.rm = T)
  xmax <- quantile(banddf[, (i+2)], UL, na.rm = T)
  x <- clamp(hardinlsat[[i]], lower = xmin, upper = xmax, values = T)
  mv('x', bname)
  bandnamec[i] <- bname
}


#step 11 : plot RGB and FCC Plots
lsathardinc<- c(aerosolxc, bluexc, greenxc, redxc, NIRxc,SWR1xc, SWR2xc)
plotRGB(lsathardinc, r= 4, g=3, b=2, stretch = 'hist')
plotRGB(hardinlsat, r= 5, g=3, b=2, stretch = 'hist')
#Clamped histogram
hist(lsathardinc)

#Step 12 :create a dataframe for aditional analysis
bandcdf <- as.data.frame(hardinlsat, xy = T, geom= 'WKT')
summary(bandcdf)
bandcdf <- na.omit(bandcdf)
colnames(bandcdf) <- c('X', 'Y', bandnames,'WKT')

#step 13: compute correlation between bands
corbands <- cor(bandcdf[,3:8], method = 'spearman')
corrplot(corbands, method = 'number', type = 'lower', diag = FALSE)

#Step 14 :Reproject data to lat-lon
crsaea <- crs(lsathardinc, proj = T) #Get crs
crs84 <- 4326 #EPSG code #Define the new crs

banddf.SP <- st_as_sf(bandcdf, coords = c('X', 'Y'), crs = crsaea)
banddf.SP$XAEA = st_coordinates(banddf.SP)[,1]
banddf.SP$YAEA = st_coordinates(banddf.SP)[,2]

bandcdf.SP <- st_transform(x = banddf.SP, crs = crs84) #reproject
bandcdf.SP$Lon = st_coordinates(bandcdf.SP)[,1]
bandcdf.SP$Lat = st_coordinates(bandcdf.SP)[,2]


#Step 15 :create a csv file
bandcdf.SP <- subset(bandcdf.SP, select = -c(WKT, geometry))
write.csv(bandcdf.SP, 'bandcdf.csv', row.names = F)

# Step 16: Perform PCA with bands only
pca <- prcomp(bandcdf[,3:8], scale = TRUE)

#Step 17: Get the first two principal components

scaled_data <- scale(bandcdf[, 3:8])
pca_result <- prcomp(scaled_data)
summary(pca_result)
pcadata <- pca$x
write.csv(pcadata[,1:3], 'pcdata.csv', row.names = F)
# Step 18: Look at the first few rows of the `rotation` matrix
head(pca$rotation)
pc1 <- pca$x[ ,2]
pc2 <- pca$x[ ,3]


#Step 19: Visualize the data points in the first two principal components
ggplot(bandcdf[, 3:6], aes(x = pc1, y = pc2, color = "darkblue")) +
  geom_point() +
  labs(title = " first two principal components",
       x = "PC1", y = "PC2")

# Step 20 :perform clustering
bpcaclust <- kmeans(pcadata[,1:3], centers=5, nstart = 10)
summary(bpcaclust)
fviz_cluster(bpcaclust, data = pcadata[,1:3],
             geom = "point",
             ellipse.type = 'convex')

bpcclust <- kmeans(pcadata[,1:3],5,10)
clusts <- bpcclust$cluster
pcadat.clust <- data.frame(bandcdf.SP,clusts)
write.csv(pcadat.clust,'bandcdfclus.csv',row.names=F)
