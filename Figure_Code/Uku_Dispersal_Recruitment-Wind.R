#making Figure 4 for the uku dispersal paper
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

Recruitment_North_Wind<-glm(Rec~Mean_V, data=Wind_Data_Annual, family="poisson")
summary(Recruitment_North_Wind)
plot(Recruitment_North_Wind)

Recruitment_Ratio_North_Wind<-glm(Rec~Mean_V, data=Wind_Data_Annual, family="poisson")
summary(Recruitment_Ratio_North_Wind)

Wind_Data_Annual$North_Wind_Rec_Predictions<-predict(Recruitment_North_Wind, Wind_Data_Annual, type="response")
cor(Wind_Data_Annual$Rec[1:31], Wind_Data_Annual$North_Wind_Rec_Predictions[1:31])^2


Rec_SSB_North_Wind<-glm(Rec_SSB_Ratio~Mean_V, data=Wind_Data_Annual, family=gaussian(link="log"))
summary(Rec_SSB_North_Wind)

Wind_Data_Annual$North_Wind_Rec_SSB_Predictions<-predict(Rec_SSB_North_Wind, Wind_Data_Annual, type="response")
cor(Wind_Data_Annual$Rec[1:31], Wind_Data_Annual$North_Wind_Rec_SSB_Predictions[1:31])^2


Wind_Data_Annual$Recruits_per_SSB<-((Wind_Data_Annual$Rec)*1000)/(Wind_Data_Annual$SSB)
hist(log(Wind_Data_Annual$Recruits_per_SSB))


Rec_SSB_Mean_V_GLM<-glm(Recruits_per_SSB~Mean_V, data=Wind_Data_Annual, family=gaussian(link="log"))
summary(Rec_SSB_Mean_V_GLM)

plot(Rec_SSB_Mean_V_GLM)
PseudoR2(Rec_SSB_Mean_V_GLM, which="all")


Rec_SSB_Mean_V_Dat<-seq(-3,-1.2, by=0.001)
Rec_SSB_Mean_V_Dat<-as.data.frame(Rec_SSB_Mean_V_Dat)
colnames(Rec_SSB_Mean_V_Dat)<-"Mean_V"
#very very similar to betareg output, should not be an issue
Rec_SSB_Mean_V_Est<-cbind(as.data.frame(predict(Rec_SSB_Mean_V_GLM, Rec_SSB_Mean_V_Dat, se.fit=TRUE, type='link')), Rec_SSB_Mean_V_Dat)
Rec_SSB_Mean_V_Est <- mutate(Rec_SSB_Mean_V_Est, lwr = exp((fit - (2.1 * se.fit))), upr = exp((fit + se.fit))) # calculating the 95% confidence interval
Rec_SSB_Mean_V_Est$Recruits_per_SSB<-exp(Rec_SSB_Mean_V_Est$fit)

png("Rec_SSB_Mean_V_Model.png", height=5, width=6, res=300, units="in")
p<-ggplot(Wind_Data_Annual, aes(x=Mean_V, y= Recruits_per_SSB)) +
  geom_point() +geom_smooth(data = Rec_SSB_Mean_V_Est, aes(ymin = lwr, ymax = upr), stat = 'identity')  + theme_classic()+ylab("Recruits per mt of SSB")+xlab("Mean Meridional Wind Strength")+theme(axis.text=element_text(size=14),axis.title=element_text(size=18,face="bold"), plot.title = element_text(size = 20, face = "bold"))
print(p)
dev.off()

Wind_Data_Annual$Recruits<-Wind_Data_Annual$Rec
Recruits_Mean_V_GLM<-glm(Recruits~Mean_V, data=Wind_Data_Annual, family=gaussian(link="log"))
summary(Recruits_Mean_V_GLM)
plot(Recruits_Mean_V_GLM)
PseudoR2(Recruits_Mean_V_GLM, which="all")


Recruits_Mean_V_Dat<-seq(-3,-1.2, by=0.001)
Recruits_Mean_V_Dat<-as.data.frame(Recruits_Mean_V_Dat)
colnames(Recruits_Mean_V_Dat)<-"Mean_V"
#very very similar to betareg output, should not be an issue
Recruits_Mean_V_Est<-cbind(as.data.frame(predict(Recruits_Mean_V_GLM, Recruits_Mean_V_Dat, se.fit=TRUE, type='link')), Recruits_Mean_V_Dat)
Recruits_Mean_V_Est <- mutate(Recruits_Mean_V_Est, lwr = exp((fit - (2.1 * se.fit))), upr = exp((fit + se.fit))) # calculating the 95% confidence interval
Recruits_Mean_V_Est$Recruits<-exp(Recruits_Mean_V_Est$fit)

png("Recruits_Mean_V_Model.png", height=5, width=6, res=300, units="in")
p<-ggplot(Wind_Data_Annual, aes(x=Mean_V, y= Recruits)) +
  geom_point() +geom_smooth(data = Recruits_Mean_V_Est, aes(ymin = lwr, ymax = upr), stat = 'identity')  + theme_classic()+ylab("Recruits (1000s)")+xlab("Mean Meridional Wind Strength")+theme(axis.text=element_text(size=14),axis.title=element_text(size=18,face="bold"), plot.title = element_text(size = 20, face = "bold"))
print(p)
dev.off()


########wind speed############

Wind_Data_Annual$Recruits_per_SSB<-((Wind_Data_Annual$Rec)*1000)/(Wind_Data_Annual$SSB)
hist(log(Wind_Data_Annual$Recruits_per_SSB))


Rec_SSB_Mean_Wind_Speed_GLM<-glm(Recruits_per_SSB~Mean_Wind_Speed, data=Wind_Data_Annual, family=gaussian(link="log"))
summary(Rec_SSB_Mean_Wind_Speed_GLM)
Wind_Data_Annual$Wind_Speed_Rec_SSB_Predictions<-predict(Rec_SSB_Mean_Wind_Speed_GLM, Wind_Data_Annual, type="response")
cor(Wind_Data_Annual$Recruits_per_SSB[1:31], Wind_Data_Annual$Wind_Speed_Rec_SSB_Predictions[1:31])^2




plot(Rec_SSB_Mean_Wind_Speed_GLM)
PseudoR2(Rec_SSB_Mean_Wind_Speed_GLM, which="all")


Rec_SSB_Mean_Wind_Speed_Dat<-seq(6.2,8, by=0.001)
Rec_SSB_Mean_Wind_Speed_Dat<-as.data.frame(Rec_SSB_Mean_Wind_Speed_Dat)
colnames(Rec_SSB_Mean_Wind_Speed_Dat)<-"Mean_Wind_Speed"
#very very similar to betareg output, should not be an issue
Rec_SSB_Mean_Wind_Speed_Est<-cbind(as.data.frame(predict(Rec_SSB_Mean_Wind_Speed_GLM, Rec_SSB_Mean_Wind_Speed_Dat, se.fit=TRUE, type='link')), Rec_SSB_Mean_Wind_Speed_Dat)
Rec_SSB_Mean_Wind_Speed_Est <- mutate(Rec_SSB_Mean_Wind_Speed_Est, lwr = exp((fit - (2.1 * se.fit))), upr = exp((fit + se.fit))) # calculating the 95% confidence interval
Rec_SSB_Mean_Wind_Speed_Est$Recruits_per_SSB<-exp(Rec_SSB_Mean_Wind_Speed_Est$fit)

png("Rec_SSB_Mean_Wind_Speed_Model.png", height=5, width=6, res=300, units="in")
p<-ggplot(Wind_Data_Annual, aes(x=Mean_Wind_Speed, y= Recruits_per_SSB)) +
  geom_point() +geom_smooth(data = Rec_SSB_Mean_Wind_Speed_Est, aes(ymin = lwr, ymax = upr), stat = 'identity', linetype = "dashed")  + theme_classic()+ylab("Recruits per mt of SSB")+xlab("Mean Wind Strength")+theme(axis.text=element_text(size=14),axis.title=element_text(size=18,face="bold"), plot.title = element_text(size = 20, face = "bold"))
print(p)
dev.off()

Wind_Data_Annual$Recruits<-Wind_Data_Annual$Rec
Recruits_Mean_Wind_Speed_GLM<-glm(Rec~Mean_Wind_Speed, data=Wind_Data_Annual, family="poisson")
summary(Recruits_Mean_Wind_Speed_GLM)

Wind_Data_Annual$Wind_Speed_Rec_Predictions<-predict(Recruits_Mean_Wind_Speed_GLM, Wind_Data_Annual, type="response")

cor(Wind_Data_Annual$Rec[1:31], Wind_Data_Annual$Wind_Speed_Rec_Predictions[1:31])^2
plot(Recruits_Mean_Wind_Speed_GLM)
PseudoR2(Recruits_Mean_Wind_Speed_GLM, which="all")


Recruits_Mean_Wind_Speed_Dat<-seq(6.2,8, by=0.001)
Recruits_Mean_Wind_Speed_Dat<-as.data.frame(Recruits_Mean_Wind_Speed_Dat)
colnames(Recruits_Mean_Wind_Speed_Dat)<-"Mean_Wind_Speed"
#very very similar to betareg output, should not be an issue
Recruits_Mean_Wind_Speed_Est<-cbind(as.data.frame(predict(Recruits_Mean_Wind_Speed_GLM, Recruits_Mean_Wind_Speed_Dat, se.fit=TRUE, type='link')), Recruits_Mean_Wind_Speed_Dat)
Recruits_Mean_Wind_Speed_Est <- mutate(Recruits_Mean_Wind_Speed_Est, lwr = exp((fit - (2.1 * se.fit))), upr = exp((fit + se.fit))) # calculating the 95% confidence interval
Recruits_Mean_Wind_Speed_Est$Recruits<-exp(Recruits_Mean_Wind_Speed_Est$fit)

png("Recruits_Mean_Wind_Speed_Model.png", height=5, width=6, res=300, units="in")
p<-ggplot(Wind_Data_Annual, aes(x=Mean_Wind_Speed, y= Recruits)) +
  geom_point() +geom_smooth(data = Recruits_Mean_Wind_Speed_Est, aes(ymin = lwr, ymax = upr), stat = 'identity')  + theme_classic()+ylab("Recruits (1000s)")+xlab("Mean Wind Strength")+theme(axis.text=element_text(size=14),axis.title=element_text(size=18,face="bold"), plot.title = element_text(size = 20, face = "bold"))
print(p)
dev.off()