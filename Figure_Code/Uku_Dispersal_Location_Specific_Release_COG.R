#making Figure 3 for the uku dispersal paper
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
library(viridis)
library(plyr)
setwd("/home/jsuca/Uku_Dispersal_Code")

#read in the individual connectivity data
Kauai_Connectivity<-readRDS( "Kauai_Connectivity_by_Year_2020_150m_PB.rds")
Hawaii_Connectivity<-readRDS("Hawaii_Connectivity_by_Year_2020_150m_PB.rds")
Oahu_Connectivity<-readRDS("Oahu_Connectivity_by_Year_2020_150m_PB.rds")
Maui_Nui_Connectivity<-readRDS("Maui_Nui_Connectivity_by_Year_2020_150m_PB.rds")
Niihau_Connectivity<-readRDS("Niihau_Connectivity_by_Year_2020_150m_PB.rds")
Kaula_Connectivity<-readRDS("Kaula_Connectivity_by_Year_2020_150m_PB.rds")

Kauai_Connectivity$Proportion_Connected<-rowMeans(Kauai_Connectivity[,c(5:17)])
Hawaii_Connectivity$Proportion_Connected<-rowMeans(Hawaii_Connectivity[,c(5:17)])
Maui_Nui_Connectivity$Proportion_Connected<-rowMeans(Maui_Nui_Connectivity[,c(5:17)])
Niihau_Connectivity$Proportion_Connected<-rowMeans(Niihau_Connectivity[,c(5:17)])
Oahu_Connectivity$Proportion_Connected<-rowMeans(Oahu_Connectivity[,c(5:17)])
Kaula_Connectivity$Proportion_Connected<-rowMeans(Kaula_Connectivity[,c(5:17)])

setwd("/home/jsuca/Uku_Dispersal_Code/Figure_Code_Data")

Hawaii_EEZ_bathy<-raster("hi_eez_extract_grid.tiff")

Ex<-c(-158, -157, 20, 22)

Hawaii_MHI<-crop(Hawaii_EEZ_bathy, Ex)
plot(Hawaii_MHI)
Hawaii_EEZ_DF<-as.data.frame(Hawaii_MHI, xy=TRUE)

colnames(Hawaii_EEZ_DF)<-c("Lon","Lat","Depth")
world <- ne_countries(scale=10,returnclass = "sf")#generate high res coastlines 

Islands<-c("Hawaii","Maui_Nui","Oahu","Kauai","Niihau", "Kaula")
for (i in 1:length(Islands)){
  Data<-get(paste0(Islands[i],"_Connectivity"))
  Data$Round_Lon<-round(Data$V2, digits=1)
  Data$Round_Lat<-round(Data$V3, digits=1)
  Spat_Summary<-ddply(Data, .(Round_Lon, Round_Lat), summarize, Proportion_Connected=mean(Proportion_Connected))
  
  Mean_Lat<-sum(Spat_Summary$Round_Lat*Spat_Summary$Proportion_Connected)/sum(Spat_Summary$Proportion_Connected)
  Mean_Lon<-sum(Spat_Summary$Round_Lon*Spat_Summary$Proportion_Connected)/sum(Spat_Summary$Proportion_Connected)
  png(paste0("Connectivity_to_",Islands[i],"_COG.png"), height=5, width=7, units="in", res=300)
  p<-ggplot()+geom_raster(data=Hawaii_EEZ_DF,aes(x=Lon, y=Lat, fill=Depth))+ geom_sf(data = world) +geom_point(data=get(paste0(Islands[i],"_Connectivity")),aes(x=V2, y=V3, colour=Proportion_Connected), size=0.2)+scale_fill_gradient(low = "gray", high = "gray3",limits=c(-7000,100))+scale_color_viridis()+geom_point(aes(x=Mean_Lon, y=Mean_Lat), shape="\u2605", size=10,colour = "red")+ coord_sf(xlim = c(-157.88, -157.22), ylim = c(20.75,21.25)) +ylab("Latitude")+xlab("Longitude")+theme_bw()
  print(p)
  rm(Mean_Lat)
  rm(Mean_Lon)
  rm(Data)
  rm(Spat_Summary)
  dev.off()
}