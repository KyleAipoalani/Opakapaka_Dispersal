#looking at retention patterns through time 

library(ncdf4)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(raster)
library(ggplot2)
library(ggnewscale)
library(sp)
library(maptools)
library(PBSmapping)
library(matrixStats)
setwd("/home/jsuca/Uku_Dispersal_Code/Figure_Code_Data")

Hawaii_EEZ_bathy<-raster("hi_eez_extract_grid.tiff")
Ex<-extent(-163, -154, 18, 25)

Hawaii_MHI<-crop(Hawaii_EEZ_bathy, Ex)
plot(Hawaii_MHI)
Hawaii_EEZ_DF<-as.data.frame(Hawaii_MHI, xy=TRUE)
colnames(Hawaii_EEZ_DF)<-c("Lon","Lat","Depth")
world <- ne_countries(scale=10,returnclass = "sf")#generate high res coastlines 
png("MHI_Plot_150m_Isobath.png", height=5, width=7, units="in", res=300)
p<-ggplot()+geom_raster(data=Hawaii_EEZ_DF,aes(x=Lon, y=Lat, fill=Depth))+ geom_sf(data = world) +geom_contour(data = Hawaii_EEZ_DF, 
                                                                                                               aes(x = Lon, y = Lat, z = Depth),
                                                                                                               breaks = -150, color = "white", linewidth = 0.8)+ scale_fill_gradient(low = "gray", high = "gray3",limits=c(-7000,100))+coord_sf(xlim = c(-160.75, -155), ylim = c(18.5,23.5)) +ylab("Latitude")+xlab("Longitude")+theme_bw()
print(p)
dev.off()

Release_locs<-read.csv("Uku_General_Habitat_PB.csv", header=FALSE)

colnames(Release_locs)<-c("Count","Lon","Lat","Pr_1")
Ex<-c(-158, -157, 20, 22)

Hawaii_MHI<-crop(Hawaii_EEZ_bathy, Ex)
plot(Hawaii_MHI)
Hawaii_EEZ_DF<-as.data.frame(Hawaii_MHI, xy=TRUE)
colnames(Hawaii_EEZ_DF)<-c("Lon","Lat","Depth")

png("MHI_Release_Plot_Isobath.png", height=5, width=7, units="in", res=300)
p<-ggplot()+geom_raster(data=Hawaii_EEZ_DF,aes(x=Lon, y=Lat, fill=Depth))+ geom_sf(data = world) + geom_sf(data = world) +geom_point(aes(x=Release_locs$Lon, y=Release_locs$Lat, colour="red"), size=0.2)+geom_contour(data = Hawaii_EEZ_DF, 
                                                                                                                                                                                                                       aes(x = Lon, y = Lat, z = Depth),
                                                                                                                                                                                                                       breaks = -150, color = "white", linewidth = 0.8)+ scale_fill_gradient(low = "gray", high = "gray3",limits=c(-7000,100))+coord_sf(xlim = c(-157.88, -157.22), ylim = c(20.75,21.25)) +ylab("Latitude")+xlab("Longitude")+theme_bw()+theme(legend.position = "none")
print(p)
dev.off()