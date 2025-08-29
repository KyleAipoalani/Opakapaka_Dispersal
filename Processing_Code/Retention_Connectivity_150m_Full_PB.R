#making plots of the diseprsal model output 


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
Release_Locations<-read.csv("Uku_General_Habitat_PB.csv", header=FALSE)
Hawaii_Connectivity_by_Year<-matrix(, nrow=nrow(Release_Locations), ncol=13)
Yrs<-c(2008:2020)
for (i in c(1:length(Yrs))){
  Hawaii_Settlement<-readRDS(paste0("Hawaii_Settlement_",Yrs[i],"_20m_150m_PB.rds"))
  Hawaii_PLD_Summed<-colSums(Hawaii_Settlement, na.rm=FALSE)
  Hawaii_PLD_Summed[Hawaii_PLD_Summed>0]<-1
  Hawaii_PLD_Date_Loc<-colSums( matrix(Hawaii_PLD_Summed, nrow=3))
  Hawaii_PLD_Loc<-rowSums( matrix(Hawaii_PLD_Date_Loc, nrow=74544))/(3*84)
  Hawaii_Connectivity_by_Year[,i]<-Hawaii_PLD_Loc
  rm(Hawaii_PLD_Loc)
  print(paste("Finished", i, "of", length(Yrs)))}

Hawaii_Connectivity<-cbind(Release_Locations, Hawaii_Connectivity_by_Year)
saveRDS(Hawaii_Connectivity, "Hawaii_Connectivity_by_Year_2020_150m_PB.rds")

gc()

Oahu_Connectivity_by_Year<-matrix(, nrow=nrow(Release_Locations), ncol=13)
Yrs<-c(2008:2020)
for (i in c(1:length(Yrs))){
  Oahu_Settlement<-readRDS(paste0("Oahu_Settlement_",Yrs[i],"_20m_150m_PB.rds"))
  Oahu_PLD_Summed<-colSums(Oahu_Settlement, na.rm=FALSE)
  Oahu_PLD_Summed[Oahu_PLD_Summed>0]<-1
  Oahu_PLD_Date_Loc<-colSums( matrix(Oahu_PLD_Summed, nrow=5))
  Oahu_PLD_Loc<-rowSums( matrix(Oahu_PLD_Date_Loc, nrow=74544))/(3*84)
  Oahu_Connectivity_by_Year[,i]<-Oahu_PLD_Loc
  rm(Oahu_PLD_Loc)
  print(paste("Finished", i, "of", length(Yrs)))
}

Oahu_Connectivity<-cbind(Release_Locations, Oahu_Connectivity_by_Year)
saveRDS(Oahu_Connectivity, "Oahu_Connectivity_by_Year_2020_150m_PB.rds")

gc()


Kauai_Connectivity_by_Year<-matrix(, nrow=nrow(Release_Locations), ncol=13)
Yrs<-c(2008:2020)
for (i in c(1:length(Yrs))){
  Kauai_Settlement<-readRDS(paste0("Kauai_Settlement_",Yrs[i],"_20m_150m_PB.rds"))
  Kauai_PLD_Summed<-colSums(Kauai_Settlement, na.rm=FALSE)
  Kauai_PLD_Summed[Kauai_PLD_Summed>0]<-1
  Kauai_PLD_Date_Loc<-colSums( matrix(Kauai_PLD_Summed, nrow=5))
  Kauai_PLD_Loc<-rowSums( matrix(Kauai_PLD_Date_Loc, nrow=74544))/(3*84)
  Kauai_Connectivity_by_Year[,i]<-Kauai_PLD_Loc
  rm(Kauai_PLD_Loc)
  print(paste("Finished", i, "of", length(Yrs)))
}

Kauai_Connectivity<-cbind(Release_Locations, Kauai_Connectivity_by_Year)
saveRDS(Kauai_Connectivity, "Kauai_Connectivity_by_Year_2020_150m_PB.rds")

gc()

Maui_Connectivity_by_Year<-matrix(, nrow=nrow(Release_Locations), ncol=13)
Yrs<-c(2008:2020)
for (i in c(1:length(Yrs))){
  Maui_Settlement<-readRDS(paste0("Maui_Nui_Settlement_",Yrs[i],"_20m_150m_PB.rds"))
  Maui_PLD_Summed<-colSums(Maui_Settlement, na.rm=FALSE)
  Maui_PLD_Summed[Maui_PLD_Summed>0]<-1
  Maui_PLD_Date_Loc<-colSums( matrix(Maui_PLD_Summed, nrow=5))
  Maui_PLD_Loc<-rowSums( matrix(Maui_PLD_Date_Loc, nrow=74544))/(3*84)
  Maui_Connectivity_by_Year[,i]<-Maui_PLD_Loc
  rm(Maui_PLD_Loc)
  print(paste("Finished", i, "of", length(Yrs)))
}

Maui_Connectivity<-cbind(Release_Locations, Maui_Connectivity_by_Year)
saveRDS(Maui_Connectivity, "Maui_Nui_Connectivity_by_Year_2020_150m_PB.rds")

Niihau_Connectivity_by_Year<-matrix(, nrow=nrow(Release_Locations), ncol=13)
Yrs<-c(2008:2020)
for (i in c(1:length(Yrs))){
  Niihau_Settlement<-readRDS(paste0("Niihau_Settlement_",Yrs[i],"_20m_150m_PB.rds"))
  Niihau_PLD_Summed<-colSums(Niihau_Settlement, na.rm=FALSE)
  Niihau_PLD_Summed[Niihau_PLD_Summed>0]<-1
  Niihau_PLD_Date_Loc<-colSums( matrix(Niihau_PLD_Summed, nrow=5))
  Niihau_PLD_Loc<-rowSums( matrix(Niihau_PLD_Date_Loc, nrow=74544))/(3*84)
  Niihau_Connectivity_by_Year[,i]<-Niihau_PLD_Loc
  rm(Niihau_PLD_Loc)
  print(paste("Finished", i, "of", length(Yrs)))
}

Niihau_Connectivity<-cbind(Release_Locations, Niihau_Connectivity_by_Year)
saveRDS(Niihau_Connectivity, "Niihau_Connectivity_by_Year_2020_150m_PB.rds")
gc()

Kaula_Connectivity_by_Year<-matrix(, nrow=nrow(Release_Locations), ncol=13)
Yrs<-c(2008:2020)
for (i in c(1:length(Yrs))){
  Kaula_Settlement<-readRDS(paste0("Kaula_Settlement_",Yrs[i],"_20m_150m_PB.rds"))
  Kaula_PLD_Summed<-colSums(Kaula_Settlement, na.rm=FALSE)
  Kaula_PLD_Summed[Kaula_PLD_Summed>0]<-1
  Kaula_PLD_Date_Loc<-colSums( matrix(Kaula_PLD_Summed, nrow=5))
  Kaula_PLD_Loc<-rowSums( matrix(Kaula_PLD_Date_Loc, nrow=74544))/(3*84)
  Kaula_Connectivity_by_Year[,i]<-Kaula_PLD_Loc
  rm(Kaula_PLD_Loc)
  print(paste("Finished", i, "of", length(Yrs)))
}

Kaula_Connectivity<-cbind(Release_Locations, Kaula_Connectivity_by_Year)
saveRDS(Kaula_Connectivity, "Kaula_Connectivity_by_Year_2020_150m_PB.rds")