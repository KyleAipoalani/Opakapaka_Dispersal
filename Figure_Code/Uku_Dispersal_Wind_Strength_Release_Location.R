#making Figure for the uku dispersal paper
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

Wind_Data_Annual<-read.csv("Yearly_Maui_Nui_Winds_Trades_High_Res_PB.csv", header=TRUE)

Islands<-c("Hawaii","Maui_Nui","Oahu","Kauai","Niihau", "Kaula")

Mean_Wind_Annual<-Wind_Data_Annual[21:33,5]
for (i in 1:length(Islands)){
  Corr_Coeff_Maui_Nui_Retention<-matrix(,nrow(get(paste0(Islands[i],"_Connectivity"))),3)
  
  for (ii in 1:nrow(Corr_Coeff_Maui_Nui_Retention)){
    Resp<-get(paste0(Islands[i],"_Connectivity"))[ii,5:17]
    Resp1<-t(Resp)
    DF<-as.data.frame(cbind(Resp1, Mean_Wind_Annual))
    names(DF)[1] <- 'Retention'
    X<-cor.test(DF$Mean_Wind_Annual, DF$Retention)
    Corr_Coeff_Maui_Nui_Retention[ii,1]<-X$estimate
    Corr_Coeff_Maui_Nui_Retention[ii,2]<-X$p.value
    Corr_Coeff_Maui_Nui_Retention[ii,3]<-mean(DF$Retention, na.rm=T)
    #print(paste("Finished",ii,"of",nrow(Corr_Coeff_Maui_Nui_Retention) ))
    
  }
  
  
  Release_Locations_Maui_Nui<-read.csv("Uku_General_Habitat_PB.csv", header=FALSE)
  Release_Locations_Maui_Nui$Corr<-Corr_Coeff_Maui_Nui_Retention[,1]
  Release_Locations_Maui_Nui$Sig<-Corr_Coeff_Maui_Nui_Retention[,2]
  Release_Locations_Maui_Nui$Mean<-Corr_Coeff_Maui_Nui_Retention[,3]
  Release_Locations_Sig_Maui_Nui<-Release_Locations_Maui_Nui[Release_Locations_Maui_Nui$Sig<=0.05 & Release_Locations_Maui_Nui$Mean>=0.01,]
  
  
  
  
  
  
  png(paste0("Correlation_with_Wind_Strength_",Islands[i],".png"), height=5, width=7, units="in", res=300)
  p<-ggplot()+geom_raster(data=Hawaii_EEZ_DF,aes(x=Lon, y=Lat, fill=Depth))+ geom_sf(data = world) +geom_point(data= Release_Locations_Sig_Maui_Nui,aes(x=V2, y=V3, colour=Corr), size=0.2)+scale_colour_gradient2(limits=c(-1,1),low ="blue",mid = "white", high = "red")+ scale_fill_gradient(low = "gray", high = "gray3",limits=c(-7000,100))+ coord_sf(xlim = c(-157.88, -157.22), ylim = c(20.75,21.25)) +ylab("Latitude")+xlab("Longitude")+theme_bw()
  print(p)
  dev.off()
}