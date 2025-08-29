#making Figure  for the uku dispersal paper
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
library(betareg)
library(mgcv)
library(cowplot)
library(DescTools)

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



Hawaii_Total_Connections_by_Year<-Hawaii_Connectivity[,c(5:17)]*3*84

Hawaii_Annual_Connections<-cbind(2008:2020,colSums(Hawaii_Total_Connections_by_Year))

Maui_Nui_Total_Connections_by_Year<-Maui_Nui_Connectivity[,c(5:17)]*3*84

Maui_Nui_Annual_Connections<-cbind(2008:2020,colSums(Maui_Nui_Total_Connections_by_Year))
plot(Maui_Nui_Annual_Connections[,1], Maui_Nui_Annual_Connections[,2])

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


Total_Connections_by_Year<-rbind(Hawaii_Annual_Connections, Maui_Nui_Annual_Connections, Oahu_Annual_Connections, Kauai_Annual_Connections, Niihau_Annual_Connections, Kaula_Annual_Connections)
colnames(Total_Connections_by_Year)<-c("Year","Connections")
Total_Connections_by_Year<-as.data.frame(Total_Connections_by_Year)
Total_Connections_Yearly<-ddply(Total_Connections_by_Year, .(Year), summarize, Total_Connections=sum(Connections))






setwd("/home/jsuca/Uku_Dispersal_Code/Figure_Code_Data")



Wind_Data_Annual<-read.csv("Yearly_Maui_Nui_Winds_Trades_High_Res_PB.csv", header=TRUE)

Wind_Data_Merged<-merge(Wind_Data_Annual, Total_Connections_Yearly, by=c("Year"))


Wind_Data_Merged$Prop_Connected<-Wind_Data_Merged$Total_Connections/18785088


hist(Wind_Data_Merged$Rec)

Connections_North_Wind<-betareg(Prop_Connected~Mean_V, data=Wind_Data_Merged)
summary(Connections_North_Wind)


Connections_Beta = gam(Prop_Connected ~ Mean_V, family=betar(link="logit"), data = Wind_Data_Merged)
summary(Connections_Beta)
Connections_Est<-cbind(as.data.frame(predict(Connections_Beta, Wind_Data_Merged, se.fit=TRUE, type='link')), Wind_Data_Merged)
Connections_Est <- mutate(Connections_Est, lwr = Connections_Beta$family$linkinv( fit - 2.1 * se.fit), upr = Connections_Beta$family$linkinv(fit + 2.1 * se.fit)) # calculating the 95% confidence interval
Connections_Est$Prop_Connected<-Connections_Beta$family$linkinv(Connections_Est$fit)


png("Betaregression_Connections_North_Winds.png", height=5, width=7, units="in", res=300)
P<-ggplot(Wind_Data_Merged, aes(x = Mean_V, y = Prop_Connected)) +
  geom_point(size = 4, shape = 19) +
  scale_fill_grey() +
  geom_smooth(data = Connections_Est, aes(ymin = lwr, ymax = upr), stat = 'identity')+
  theme_cowplot()+theme(legend.position = "none")+ylab("Proportion Connected")+xlab("Mean Meridional Wind Strength")
print(P)
dev.off()


########wind speed############

Connections_Wind_Speed<-betareg(Prop_Connected~Mean_Wind_Speed, data=Wind_Data_Merged)
summary(Connections_Wind_Speed)


Connections_Beta = gam(Prop_Connected ~ Mean_Wind_Speed, family=betar(link="logit"), data = Wind_Data_Merged)
summary(Connections_Beta)
Connections_Est<-cbind(as.data.frame(predict(Connections_Beta, Wind_Data_Merged, se.fit=TRUE, type='link')), Wind_Data_Merged)
Connections_Est <- mutate(Connections_Est, lwr = Connections_Beta$family$linkinv( fit - 2.1 * se.fit), upr = Connections_Beta$family$linkinv(fit + 2.1 * se.fit)) # calculating the 95% confidence interval
Connections_Est$Prop_Connected<-Connections_Beta$family$linkinv(Connections_Est$fit)


png("Betaregression_Connections_Mean_Wind_Speed.png", height=5, width=7, units="in", res=300)
P<-ggplot(Wind_Data_Merged, aes(x = Mean_Wind_Speed, y = Prop_Connected)) +
  geom_point(size = 4, shape = 19) +
  scale_fill_grey() +
  geom_smooth(data = Connections_Est, aes(ymin = lwr, ymax = upr), stat = 'identity')+
  theme_cowplot()+theme(legend.position = "none")+ylab("Proportion Connected")+xlab("Mean Wind Strength")
print(P)
dev.off()