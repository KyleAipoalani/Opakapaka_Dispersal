
#premilinary look at connectivity across the first 50 days


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


setwd("/home/jsuca/Uku_Dispersal_Code/Scripts/2019/20m/Penguin_Bank")
Output1<-nc_open("PB_Uku_20m_LL_2019.nc")

Lat<-ncvar_get(Output1, "lat")
Lon<-ncvar_get(Output1, "lon")                                                    

Lat_Settlment_Window<-Lat[52:94,c(seq(1, ncol(Lat), by=3))]
Lon_Settlment_Window<-Lon[52:94,c(seq(1, ncol(Lon), by=3))]

rm(Lat)
rm(Lon)
gc()

Pts<-as.data.frame(cbind(Lon=as.vector(Lon_Settlment_Window),Lat= as.vector(Lat_Settlment_Window)))
Pts$Count<-1
Pts$Lat<-round(Pts$Lat, digits=2)
Pts$Lon<-round(Pts$Lon, digits=2)

library(plyr)
library(ggplot2)
Binned_Estimates<-ddply(Pts, .(Lon, Lat), summarise,Particles=sum(Count))

Binned_Estimates$Particle_Density<-log10(Binned_Estimates$Particles)
setwd("/home/jsuca/Uku_Dispersal_Code/Figure_Code_Data")
world <- ne_countries(scale=10,returnclass = "sf")#generate high res coastlines 


png("Particle_Density_Settlement_Stages_2019_Maui_Nui.png", height=6, width=8, units="in", res=300)
ggplot()+coord_fixed(ratio = 1)+geom_raster(data= Binned_Estimates,aes(x=Lon, y=Lat, fill=Particle_Density))+ geom_sf(data = world)+scale_fill_viridis_c( guide = guide_colourbar(title="Particle Density"))+
  theme_bw()+geom_sf()+ coord_sf(xlim = c(-158.7, -156), ylim = c(20,22))+
  theme(legend.title=element_text(size=16),legend.text=element_text(size=14),legend.direction = "vertical", legend.box = "vertical")+ggtitle("Particle Density 2019")
dev.off()