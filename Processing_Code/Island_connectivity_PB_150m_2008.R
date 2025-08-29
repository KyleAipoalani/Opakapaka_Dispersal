
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


setwd("/home/jsuca/Uku_Dispersal_Code/Scripts/2008/20m/Penguin_Bank")
Output1<-nc_open("PB_Uku_20m_LL_2008.nc")

Lat<-ncvar_get(Output1, "lat")
Lon<-ncvar_get(Output1, "lon")                                                    

Lat_Settlment_Window<-Lat[52:94,]
Lon_Settlment_Window<-Lon[52:94,]

Lat_Start<-Lat[1,]
Lon_Start<-Lon[1,]
rm(Lat)
rm(Lon)
gc()
setwd("/home/jsuca/Uku_Dispersal_Code")

MHI_Buffer<-readRDS("MHI_Depth_150m_Buffer.rds")


#plot(MHI_Buffer)
MHI_Buffer[is.na(MHI_Buffer)]<-0


#now look into the connectivity among locations 

MHI_Polygons<-rasterToContour(MHI_Buffer)
var_temp <- SpatialLines2PolySet(MHI_Polygons)
sp_try <- PolySet2SpatialPolygons(var_temp)
Hawaii_Coast_polygons<- sp_try[1, 1] # I took the first polygon from my Shapefile




Island_names<-c("Kaula", "Niihau", "Kauai", "Oahu", "Maui_Nui","ōpala1","ōpalu2", "Hawaii")

for ( i in 1: length(Island_names)){
  assign(paste0(Island_names[i],"_150m_Coords"), Hawaii_Coast_polygons@polygons[[1]]@Polygons[[i]]@coords)
}

Niihau_Settlement<-matrix(, nrow=nrow(Lat_Settlment_Window), ncol=ncol(Lat_Settlment_Window))

for ( i in 1: nrow(Lat_Settlment_Window)){
  Niihau_Settlement[i,]<-point.in.polygon(Lon_Settlment_Window[i,], Lat_Settlment_Window[i,], Niihau_150m_Coords[,1],Niihau_150m_Coords[,2] )
}
saveRDS(Niihau_Settlement, "Niihau_Settlement_2008_20m_150m_PB.rds")
gc()



Hawaii_Settlement<-matrix(, nrow=nrow(Lat_Settlment_Window), ncol=ncol(Lat_Settlment_Window))

for ( i in 1: nrow(Lat_Settlment_Window)){
  Hawaii_Settlement[i,]<-point.in.polygon(Lon_Settlment_Window[i,], Lat_Settlment_Window[i,], Hawaii_150m_Coords[,1],Hawaii_150m_Coords[,2] )
}
saveRDS(Hawaii_Settlement, "Hawaii_Settlement_2008_20m_150m_PB.rds")

gc()
#need to look into parsing these now by release location to make a plot that is interesting and relevant
#save.image("Preliminary_Settlement_Exploration_2008_20m.RData", safe=FALSE)

Oahu_Settlement<-matrix(, nrow=nrow(Lat_Settlment_Window), ncol=ncol(Lat_Settlment_Window))

for ( i in 1: nrow(Lat_Settlment_Window)){
  Oahu_Settlement[i,]<-point.in.polygon(Lon_Settlment_Window[i,], Lat_Settlment_Window[i,], Oahu_150m_Coords[,1],Oahu_150m_Coords[,2] )
}
saveRDS(Oahu_Settlement, "Oahu_Settlement_2008_20m_150m_PB.rds")

gc()


Kauai_Settlement<-matrix(, nrow=nrow(Lat_Settlment_Window), ncol=ncol(Lat_Settlment_Window))

for ( i in 1: nrow(Lat_Settlment_Window)){
  Kauai_Settlement[i,]<-point.in.polygon(Lon_Settlment_Window[i,], Lat_Settlment_Window[i,], Kauai_150m_Coords[,1],Kauai_150m_Coords[,2] )
}
saveRDS(Kauai_Settlement, "Kauai_Settlement_2008_20m_150m_PB.rds")

gc()




Maui_Nui_Settlement<-matrix(, nrow=nrow(Lat_Settlment_Window), ncol=ncol(Lat_Settlment_Window))

for ( i in 1: nrow(Lat_Settlment_Window)){
  Maui_Nui_Settlement[i,]<-point.in.polygon(Lon_Settlment_Window[i,], Lat_Settlment_Window[i,], Maui_Nui_150m_Coords[,1],Maui_Nui_150m_Coords[,2] )
}
saveRDS(Maui_Nui_Settlement, "Maui_Nui_Settlement_2008_20m_150m_PB.rds")

gc()


Kaula_Settlement<-matrix(, nrow=nrow(Lat_Settlment_Window), ncol=ncol(Lat_Settlment_Window))

for ( i in 1: nrow(Lat_Settlment_Window)){
  Kaula_Settlement[i,]<-point.in.polygon(Lon_Settlment_Window[i,], Lat_Settlment_Window[i,], Kaula_150m_Coords[,1],Kaula_150m_Coords[,2] )
}
saveRDS(Kaula_Settlement, "Kaula_Settlement_2008_20m_150m_PB.rds")

gc()

#now generate a connectivity matrix 