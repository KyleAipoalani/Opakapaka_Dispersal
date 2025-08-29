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

#read in the individual connectivity data
Kauai_Connectivity<-readRDS( "Kauai_Connectivity_by_Year_2020_150m_PB.rds")
Hawaii_Connectivity<-readRDS("Hawaii_Connectivity_by_Year_2020_150m_PB.rds")
Oahu_Connectivity<-readRDS("Oahu_Connectivity_by_Year_2020_150m_PB.rds")
Maui_Connectivity<-readRDS("Maui_Nui_Connectivity_by_Year_2020_150m_PB.rds")
Niihau_Connectivity<-readRDS("Niihau_Connectivity_by_Year_2020_150m_PB.rds")
Kaula_Connectivity<-readRDS("Kaula_Connectivity_by_Year_2020_150m_PB.rds")
#make plot of the connectivity by location to each major domain





Kauai_Connectivity$Total_Connected<-(rowSums(Kauai_Connectivity[,c(5:17)])*(3*84))
Hawaii_Connectivity$Total_Connected<-(rowSums(Hawaii_Connectivity[,c(5:17)])*(3*84))
Maui_Connectivity$Total_Connected<-(rowSums(Maui_Connectivity[,c(5:17)])*(3*84))
Niihau_Connectivity$Total_Connected<-(rowSums(Niihau_Connectivity[,c(5:17)])*(3*84))
Oahu_Connectivity$Total_Connected<-(rowSums(Oahu_Connectivity[,c(5:17)])*(3*84))
Kaula_Connectivity$Total_Connected<-(rowSums(Kaula_Connectivity[,c(5:17)])*(3*84))





Hawaii_Total_Connections_by_Year<-Hawaii_Connectivity[,c(5:17)]*3*84

Hawaii_Annual_Connections<-cbind(2008:2020,colSums(Hawaii_Total_Connections_by_Year))

Maui_Total_Connections_by_Year<-Maui_Connectivity[,c(5:17)]*3*84

Maui_Annual_Connections<-cbind(2008:2020,colSums(Maui_Total_Connections_by_Year))
plot(Maui_Annual_Connections[,1], Maui_Annual_Connections[,2])

Oahu_Total_Connections_by_Year<-Oahu_Connectivity[,c(5:17)]*3*84
Oahu_Annual_Connections<-cbind(2008:2020,colSums(Oahu_Total_Connections_by_Year))
plot(Oahu_Annual_Connections[,1], Oahu_Annual_Connections[,2])
Kauai_Total_Connections_by_Year<-Kauai_Connectivity[,c(5:17)]*3*84
Kauai_Annual_Connections<-cbind(2008:2020,colSums(Kauai_Total_Connections_by_Year))
plot(Kauai_Annual_Connections[,1], Kauai_Annual_Connections[,2])

Niihau_Total_Connections_by_Year<-Niihau_Connectivity[,c(5:17)]*3*84
Niihau_Annual_Connections<-cbind(2008:2020,colSums(Niihau_Total_Connections_by_Year))
plot(Niihau_Annual_Connections[,1], Niihau_Annual_Connections[,2]/max(Niihau_Annual_Connections[,2]))

Kaula_Total_Connections_by_Year<-Kaula_Connectivity[,c(5:17)]*3*84
Kaula_Annual_Connections<-cbind(2008:2020,colSums(Kaula_Total_Connections_by_Year))
plot(Kaula_Annual_Connections[,1], Kaula_Annual_Connections[,2]/max(Kaula_Annual_Connections[,2]))

Total_Connections_by_Year<-rbind(Hawaii_Annual_Connections, Maui_Annual_Connections, Oahu_Annual_Connections, Kauai_Annual_Connections, Niihau_Annual_Connections)
colnames(Total_Connections_by_Year)<-c("Year","Connections")
Total_Connections_by_Year<-as.data.frame(Total_Connections_by_Year)
Total_Connections_Yearly<-ddply(Total_Connections_by_Year, .(Year), summarize, Total_Connections=sum(Connections))




png("Connections_by_Island_Trend_No_MHI_PB.png", height=5, width=8, units="in", res=300)
plot(2008:2020, Hawaii_Annual_Connections[,2]/max(Hawaii_Annual_Connections[,2]), type="l", lwd=2, ylim=c(0,1), xlab=c("Year"), ylab=c("Std. Received"))
par(new=T)
plot(2008:2020, Maui_Annual_Connections[,2]/max(Maui_Annual_Connections[,2]), type="l", lwd=2, ylim=c(0,1),col="red", xlab=c("Year"), ylab=c("Std. Received"))
par(new=T)
plot(2008:2020, Oahu_Annual_Connections[,2]/max(Oahu_Annual_Connections[,2]), type="l", lwd=2, ylim=c(0,1),col="blue", xlab=c("Year"), ylab=c("Std. Received"))
par(new=T)
plot(2008:2020, Kauai_Annual_Connections[,2]/max(Kauai_Annual_Connections[,2]), type="l", lwd=2, ylim=c(0,1),col="green", xlab=c("Year"), ylab=c("Std. Received"))
par(new=T)
plot(2008:2020,  Niihau_Annual_Connections[,2]/max(Niihau_Annual_Connections[,2]), type="l",lty=2, ylim=c(0,1),col="brown", xlab=c("Year"), ylab=c("Std. Received"))
par(new=T)
plot(2008:2020,  Kaula_Annual_Connections[,2]/max(Kaula_Annual_Connections[,2]), type="l",lty=3, ylim=c(0,1), xlab=c("Year"), ylab=c("Std. Received"))
#par(new=T)
#plot(2008:2020,  Total_Connections_Yearly$Total_Connections/max(Total_Connections_Yearly$Total_Connections), type="l",lwd=6, ylim=c(0.1,1), xlab=c("Year"), ylab=c("Std. Received"))
dev.off()

png("Connections_by_Island_Trend_Zoomed_PB.png", height=5, width=8, units="in", res=300)
plot(2008:2020, Hawaii_Annual_Connections[,2]/max(Hawaii_Annual_Connections[,2]), type="l", lwd=2, ylim=c(0,1), xlab=c("Year"), ylab=c("Std. Received"))
par(new=T)
plot(2008:2020, Maui_Annual_Connections[,2]/max(Maui_Annual_Connections[,2]), type="l", lwd=2, ylim=c(0,1),col="red", xlab=c("Year"), ylab=c("Std. Received"))
par(new=T)
plot(2008:2020, Oahu_Annual_Connections[,2]/max(Oahu_Annual_Connections[,2]), type="l", lwd=2, ylim=c(0,1),col="blue", xlab=c("Year"), ylab=c("Std. Received"))
par(new=T)
plot(2008:2020, Kauai_Annual_Connections[,2]/max(Kauai_Annual_Connections[,2]), type="l", lwd=2, ylim=c(0,1),col="green", xlab=c("Year"), ylab=c("Std. Received"))
par(new=T)
plot(2008:2020,  Niihau_Annual_Connections[,2]/max(Niihau_Annual_Connections[,2]), type="l",lty=2, ylim=c(0,1),col="brown", xlab=c("Year"), ylab=c("Std. Received"))
par(new=T)
plot(2008:2020,  Total_Connections_Yearly$Total_Connections/max(Total_Connections_Yearly$Total_Connections), type="l",lwd=6, ylim=c(0,1), xlab=c("Year"), ylab=c("Std. Received"))
par(new=T)
plot(2008:2020,  Kaula_Annual_Connections[,2]/max(Kaula_Annual_Connections[,2]), type="l",lty=3, ylim=c(0,1), xlab=c("Year"), ylab=c("Std. Received"))

dev.off()


png("Connections_by_Island_Trend_Zoomed_PB_Prop_Releases.png", height=5, width=8, units="in", res=300)
plot(2008:2020, Hawaii_Annual_Connections[,2]/18785088, type="l", lwd=2, ylim=c(0,0.4), xlab=c("Year"), ylab=c("Prop. Connected"))
par(new=T)
plot(2008:2020, Maui_Annual_Connections[,2]/18785088, type="l", lwd=2, ylim=c(0,0.4),col="red", xlab=c("Year"), ylab=c("Prop. Connected"))
par(new=T)
plot(2008:2020, Oahu_Annual_Connections[,2]/18785088, type="l", lwd=2, ylim=c(0,0.4),col="blue", xlab=c("Year"), ylab=c("Prop. Connected"))
par(new=T)
plot(2008:2020, Kauai_Annual_Connections[,2]/18785088, type="l", lwd=2, ylim=c(0,0.4),col="green", xlab=c("Year"), ylab=c("Prop. Connected"))
par(new=T)
plot(2008:2020,  Niihau_Annual_Connections[,2]/18785088, type="l",lty=2, ylim=c(0,0.4),col="brown", xlab=c("Year"), ylab=c("Prop. Connected"))
par(new=T)
plot(2008:2020,  Kaula_Annual_Connections[,2]/18785088, type="l",lty=3, ylim=c(0,0.4), xlab=c("Year"), ylab=c("Prop. Connected"))
par(new=T)
plot(2008:2020,  Total_Connections_Yearly$Total_Connections/18785088, type="l",lwd=6, ylim=c(0,0.4), xlab=c("Year"), ylab=c("Prop. Connected"))
dev.off()