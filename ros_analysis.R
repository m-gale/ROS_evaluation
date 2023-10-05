
"""
Author: Matthew Gale matthew.gale@anu.edu.au

Code for rate of spread model evaluation, including plot generation and RF modelling.
Input: df2 - filtered csv dataframe with environmental parameters and rates of spread.

"""

library(party)
library(randomForest)
library(caTools)
library(tidyverse)
library(lme4)
library(jtools)
library(dplyr)
library(missForest)
library(Boruta)
library(rfUtilities)
library(plotmo)
library(mgcv)
library(itsadug)
library(pdp)

#####

#classify veg by code and convert to factor
df2$mvg[df2$mvg2=='O']='Other'
df2$mvg[df2$mvg2=='DS']='Dry Sclerophyll'
df2$mvg[df2$mvg2=='WS']='Wet Sclerophyll'
df2$mvg[df2$mvg2=='WL']='Eucalypt Woodlands'
df2$mvg<-as.factor(df2$mvg)
df2$mvg2<-as.factor(df2$mvg2)

#windspeed unit conversion
df2$ws_net<-df2$ws_net/1000*60*60

#calc residuals
df2$res_vesta2<-df2$ros_vesta2_afdrs-df2$rate #negative is underprediction. 
df2$res_vesta1<-df2$ros_vesta_afdrs-df2$rate
df2$res_10wind<-df2$ros_10wind-df2$rate
df2$res_mk5<-df2$ros_mk5-df2$rate

#####

#summary stats

sum(df2$mvg=='Dry Sclerophyll', na.rm=TRUE)
sum(df2$mvg=='Wet Sclerophyll', na.rm=TRUE)
sum(df2$mvg=='Eucalypt Woodlands', na.rm=TRUE)
sum(df2$mvg=='Other', na.rm=TRUE)

sum(df2$mvg=='Dry Sclerophyll')+sum(df2$mvg=='Wet Sclerophyll')+sum(df2$mvg=='Eucalypt Woodlands')+sum(df2$mvg=='Other')

mean(df2$rate[df2$mvg=='Dry Sclerophyll'])
mean(df2$rate[df2$mvg=='Wet Sclerophyll'])
mean(df2$rate[df2$mvg=='Eucalypt Woodlands'])
mean(df2$rate[df2$mvg=='Other'])
mean(df2$rate)

std.error(df2$rate[df2$mvg=='Dry Sclerophyll'])
std.error(df2$rate[df2$mvg=='Wet Sclerophyll'])
std.error(df2$rate[df2$mvg=='Eucalypt Woodlands'])
std.error(df2$rate[df2$mvg=='Other'])
std.error(df2$rate)

min(df2$rate[df2$mvg=='Dry Sclerophyll'])
min(df2$rate[df2$mvg=='Wet Sclerophyll'])
min(df2$rate[df2$mvg=='Eucalypt Woodlands'])
min(df2$rate[df2$mvg=='Other'])
min(df2$rate)

max(df2$rate[df2$mvg=='Dry Sclerophyll'])
max(df2$rate[df2$mvg=='Wet Sclerophyll'])
max(df2$rate[df2$mvg=='Eucalypt Woodlands'])
max(df2$rate[df2$mvg=='Other'])
max(df2$rate)

#fire run duration
mean(df2$sincelast[df2$mvg=='Dry Sclerophyll'])
mean(df2$sincelast[df2$mvg=='Wet Sclerophyll'])
mean(df2$sincelast[df2$mvg=='Eucalypt Woodlands'])
mean(df2$sincelast[df2$mvg=='Other'])
mean(df2$sincelast)

#####

#vegetation-wise model error stats

#mae
mean(abs(df2$ros_vesta_afdrs-df2$rate)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)
mean(abs(df2$ros_vesta2_afdrs-df2$rate)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)
mean(abs(df2$ros_10wind-df2$rate)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)
mean(abs(df2$ros_mk5-df2$rate)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)

mean(abs(df2$ros_vesta_afdrs-df2$rate)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)
mean(abs(df2$ros_vesta2_afdrs-df2$rate)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)
mean(abs(df2$ros_10wind-df2$rate)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)
mean(abs(df2$ros_mk5-df2$rate)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)

mean(abs(df2$ros_vesta_afdrs-df2$rate)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)
mean(abs(df2$ros_vesta2_afdrs-df2$rate)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)
mean(abs(df2$ros_10wind-df2$rate)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)
mean(abs(df2$ros_mk5-df2$rate)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)

#mbe
mean((df2$ros_vesta_afdrs-df2$rate)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)
mean((df2$ros_vesta2_afdrs-df2$rate)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)
mean((df2$ros_10wind-df2$rate)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)
mean((df2$ros_mk5-df2$rate)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)

mean((df2$ros_vesta_afdrs-df2$rate)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)
mean((df2$ros_vesta2_afdrs-df2$rate)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)
mean((df2$ros_10wind-df2$rate)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)
mean((df2$ros_mk5-df2$rate)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)

mean((df2$ros_vesta_afdrs-df2$rate)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)
mean((df2$ros_vesta2_afdrs-df2$rate)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)
mean((df2$ros_10wind-df2$rate)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)
mean((df2$ros_mk5-df2$rate)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)

#mbpe
mean(((df2$ros_vesta_afdrs-df2$rate)/df2$rate)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)*100
mean(((df2$ros_vesta2_afdrs-df2$rate)/df2$rate)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)*100
mean(((df2$ros_10wind-df2$rate)/df2$rate)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)*100
mean(((df2$ros_mk5-df2$rate)/df2$rate)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)*100

mean(((df2$ros_vesta_afdrs-df2$rate)/df2$rate)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)*100
mean(((df2$ros_vesta2_afdrs-df2$rate)/df2$rate)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)*100
mean(((df2$ros_10wind-df2$rate)/df2$rate)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)*100
mean(((df2$ros_mk5-df2$rate)/df2$rate)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)*100

mean(((df2$ros_vesta_afdrs-df2$rate)/df2$rate)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)*100
mean(((df2$ros_vesta2_afdrs-df2$rate)/df2$rate)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)*100
mean(((df2$ros_10wind-df2$rate)/df2$rate)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)*100
mean(((df2$ros_mk5-df2$rate)/df2$rate)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)*100

#rmse
sqrt(mean(((df2$rate - df2$ros_vesta_afdrs)^2)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE))
sqrt(mean(((df2$rate - df2$ros_vesta2_afdrs)^2)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE))
sqrt(mean(((df2$rate - df2$ros_10wind)^2)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE))
sqrt(mean(((df2$rate - df2$ros_mk5)^2)[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE))

sqrt(mean(((df2$rate - df2$ros_vesta_afdrs)^2)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE))
sqrt(mean(((df2$rate - df2$ros_vesta2_afdrs)^2)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE))
sqrt(mean(((df2$rate - df2$ros_10wind)^2)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE))
sqrt(mean(((df2$rate - df2$ros_mk5)^2)[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE))

sqrt(mean(((df2$rate - df2$ros_vesta_afdrs)^2)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE))
sqrt(mean(((df2$rate - df2$ros_vesta2_afdrs)^2)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE))
sqrt(mean(((df2$rate - df2$ros_10wind)^2)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE))
sqrt(mean(((df2$rate - df2$ros_mk5)^2)[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE))

#35%
sum((df2$rate*1.35>df2$ros_vesta_afdrs) & (df2$rate*0.65<df2$ros_vesta_afdrs)==TRUE & 
      (df2$mvg=='Dry Sclerophyll'), na.rm=TRUE)/sum((df2$ros_vesta_afdrs>=0) &  
                                                      (df2$mvg=='Dry Sclerophyll'), na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_vesta2_afdrs) & (df2$rate*0.65<df2$ros_vesta2_afdrs)==TRUE & 
      (df2$mvg=='Dry Sclerophyll'), na.rm=TRUE)/sum((df2$ros_vesta2_afdrs>=0) &  
                                                      (df2$mvg=='Dry Sclerophyll'), na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_10wind) & (df2$rate*0.65<df2$ros_10wind)==TRUE & 
      (df2$mvg=='Dry Sclerophyll'), na.rm=TRUE)/sum((df2$ros_10wind>=0) &  
                                                      (df2$mvg=='Dry Sclerophyll'), na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_mk5) & (df2$rate*0.65<df2$ros_mk5)==TRUE & 
      (df2$mvg=='Dry Sclerophyll'), na.rm=TRUE)/sum((df2$ros_mk5>=0) &  
                                                      (df2$mvg=='Dry Sclerophyll'), na.rm=TRUE)*100

sum((df2$rate*1.35>df2$ros_vesta_afdrs) & (df2$rate*0.65<df2$ros_vesta_afdrs)==TRUE & 
      (df2$mvg=='Wet Sclerophyll'), na.rm=TRUE)/sum((df2$ros_vesta_afdrs>=0) &  
                                                      (df2$mvg=='Wet Sclerophyll'), na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_vesta2_afdrs) & (df2$rate*0.65<df2$ros_vesta2_afdrs)==TRUE & 
      (df2$mvg=='Wet Sclerophyll'), na.rm=TRUE)/sum((df2$ros_vesta2_afdrs>=0) &  
                                                      (df2$mvg=='Wet Sclerophyll'), na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_10wind) & (df2$rate*0.65<df2$ros_10wind)==TRUE & 
      (df2$mvg=='Wet Sclerophyll'), na.rm=TRUE)/sum((df2$ros_10wind>=0) &  
                                                      (df2$mvg=='Wet Sclerophyll'), na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_mk5) & (df2$rate*0.65<df2$ros_mk5)==TRUE & 
      (df2$mvg=='Wet Sclerophyll'), na.rm=TRUE)/sum((df2$ros_mk5>=0) &  
                                                      (df2$mvg=='Wet Sclerophyll'), na.rm=TRUE)*100


sum((df2$rate*1.35>df2$ros_vesta_afdrs) & (df2$rate*0.65<df2$ros_vesta_afdrs)==TRUE & 
      (df2$mvg=='Eucalypt Woodlands'), na.rm=TRUE)/sum((df2$ros_vesta_afdrs>=0) &  
                                                      (df2$mvg=='Eucalypt Woodlands'), na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_vesta2_afdrs) & (df2$rate*0.65<df2$ros_vesta2_afdrs)==TRUE & 
      (df2$mvg=='Eucalypt Woodlands'), na.rm=TRUE)/sum((df2$ros_vesta2_afdrs>=0) &  
                                                      (df2$mvg=='Eucalypt Woodlands'), na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_10wind) & (df2$rate*0.65<df2$ros_10wind)==TRUE & 
      (df2$mvg=='Eucalypt Woodlands'), na.rm=TRUE)/sum((df2$ros_10wind>=0) &  
                                                      (df2$mvg=='Eucalypt Woodlands'), na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_mk5) & (df2$rate*0.65<df2$ros_mk5)==TRUE & 
      (df2$mvg=='Eucalypt Woodlands'), na.rm=TRUE)/sum((df2$ros_mk5>=0) &  
                                                      (df2$mvg=='Eucalypt Woodlands'), na.rm=TRUE)*100

mean(df2$temp[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)
mean(df2$temp[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)
mean(df2$temp[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)

mean(df2$rh[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)
mean(df2$rh[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)
mean(df2$rh[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)

mean(df2$ws_net[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)
mean(df2$ws_net[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)
mean(df2$ws_net[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)

mean(df2$st_dist[df2$mvg=='Dry Sclerophyll'], na.rm=TRUE)
mean(df2$st_dist[df2$mvg=='Wet Sclerophyll'], na.rm=TRUE)
mean(df2$st_dist[df2$mvg=='Eucalypt Woodlands'], na.rm=TRUE)


#####

#overall stats across veg types for ERA-5 weather and BoM weather ROS predictions

#MAE

mean(abs(df2$ros_vesta2_afdrs-df2$rate))
mean(abs(df2$ros_vesta_afdrs-df2$rate), na.rm=TRUE)
mean(abs(df2$ros_vesta_afdrs_era-df2$rate), na.rm=TRUE)
mean(abs(df2$ros_vesta2_afdrs_era-df2$rate))

mean(abs(df2$ros_10wind-df2$rate))
mean(abs(df2$ros_10wind_era-df2$rate))

mean(abs(df2$ros_mk5-df2$rate), na.rm=TRUE)

#MBE

mean(df2$ros_vesta2_afdrs-df2$rate)
mean(df2$ros_vesta_afdrs-df2$rate, na.rm=TRUE)
mean(df2$ros_vesta_afdrs_era-df2$rate, na.rm=TRUE)
mean(df2$ros_vesta2_afdrs_era-df2$rate)

mean(df2$ros_10wind-df2$rate)
mean(df2$ros_10wind_era-df2$rate)

mean((df2$ros_mk5-df2$rate),na.rm=TRUE)

#MAPE

mean(abs((df2$ros_vesta2_afdrs-df2$rate)/df2$rate))*100
mean(abs((df2$ros_vesta_afdrs-df2$rate)/df2$rate))*100
mean(abs((df2$ros_vesta2_afdrs_era-df2$rate)/df2$rate))*100
mean(abs((df2$ros_vesta_afdrs_era-df2$rate)/df2$rate))*100

mean(abs((df2$ros_10wind-df2$rate)/df2$rate))*100
mean(abs((df2$ros_10wind_era-df2$rate)/df2$rate))*100

mean(abs(df2$rate-df2$ros_mk5)/df2$rate*100, na.rm=TRUE)

#MBPE

mean((df2$ros_vesta2_afdrs-df2$rate)/df2$rate)*100
mean((df2$ros_vesta_afdrs-df2$rate)/df2$rate, na.rm=TRUE)*100

mean((df2$ros_vesta2_afdrs_era-df2$rate)/df2$rate)*100
mean((df2$ros_vesta_afdrs_era-df2$rate)/df2$rate, na.rm=TRUE)*100

mean((df2$ros_10wind-df2$rate)/df2$rate)*100
mean((df2$ros_10wind_era-df2$rate)/df2$rate)*100

mean((df2$ros_mk5-df2$rate)/df2$rate*100, na.rm=TRUE)

#some bias factors reported in manuscript 
mean((df2$ros_mk5/df2$rate)[df2$rate>3000],na.rm=TRUE)
mean((df2$ros_mk5/df2$rate),na.rm=TRUE)
mean((df2$ros_10wind/df2$rate),na.rm=TRUE)
mean((df2$ros_vesta2_afdrs/df2$rate),na.rm=TRUE)
mean((df2$ros_vesta_afdrs/df2$rate),na.rm=TRUE)

#RMSE

sqrt(mean((df2$rate - df2$ros_vesta2_afdrs)^2))
sqrt(mean((df2$rate - df2$ros_vesta_afdrs)^2, na.rm=TRUE))
sqrt(mean((df2$rate - df2$ros_vesta2_afdrs_era)^2))
sqrt(mean((df2$rate - df2$ros_vesta_afdrs_era)^2, na.rm=TRUE))

sqrt(mean((df2$rate - df2$ros_10wind)^2))
sqrt(mean((df2$rate - df2$ros_10wind_era)^2))

sqrt(mean((df2$rate - df2$ros_mk5)^2, na.rm=TRUE))

#%in 35%

sum((df2$rate*1.35>df2$ros_vesta_afdrs) & (df2$rate*0.65<df2$ros_vesta_afdrs)==TRUE, na.rm=TRUE)/sum(df2$ros_vesta_afdrs>=0, na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_vesta2_afdrs) & (df2$rate*0.65<df2$ros_vesta2_afdrs)==TRUE, na.rm=TRUE)/sum(df2$ros_vesta2_afdrs>=0, na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_10wind) & (df2$rate*0.65<df2$ros_10wind)==TRUE, na.rm=TRUE)/sum(df2$ros_10wind>=0, na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_mk5) & (df2$rate*0.65<df2$ros_mk5)==TRUE, na.rm=TRUE)/sum(df2$ros_mk>=0, na.rm=TRUE)*100

sum((df2$rate*1.35>df2$ros_vesta_afdrs_era) & (df2$rate*0.65<df2$ros_vesta_afdrs_era)==TRUE, na.rm=TRUE)/sum(df2$ros_vesta_afdrs_era>=0, na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_vesta2_afdrs_era) & (df2$rate*0.65<df2$ros_vesta2_afdrs_era)==TRUE, na.rm=TRUE)/sum(df2$ros_vesta2_afdrs_era>=0, na.rm=TRUE)*100
sum((df2$rate*1.35>df2$ros_10wind_era) & (df2$rate*0.65<df2$ros_10wind_era)==TRUE, na.rm=TRUE)/sum(df2$ros_10wind_era>=0, na.rm=TRUE)*100


#####

#scatter plots

loa=data.frame()

#vesta2
p1<-ggplot(df2, aes(x=rate/1000, y=ros_vesta2_afdrs/1000, colour=mvg))+
  geom_line(data=loa, aes(x=seq(0, 40, 0.01), y=seq(0, 40, 0.01)), size=0.8, linetype='solid', color='black')+
  geom_line(data=loa, aes(x=seq(0, 90, 0.01)*1.35, y=seq(0, 90, 0.01)), size=0.8, linetype='longdash', color='black')+
  geom_line(data=loa, aes(x=seq(0, 90, 0.01), y=seq(0, 90, 0.01)*1.35), size=0.8, linetype='longdash', color='black')+
  geom_point(size=4)+
  theme_bw()+
  theme(panel.border=element_blank())+
  theme(panel.grid=element_blank())+
  theme(panel.border=element_rect(colour='black', fill='transparent'))+
  theme(panel.grid=element_blank())+
  scale_y_continuous(limits=c(0, 11), breaks = seq(0, 11, by = 2), expand = c(0, 0))+
  scale_x_continuous(limits=c(0, 11), breaks = seq(0, 11, by = 2), expand = c(0, 0))+
  labs(x=expression(Observed~rate~of~spread~(km~hr^-1)), 
       y=expression(Predicted~rate~of~spread~(km~hr^-1)))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=14))+
  theme(axis.text.y= element_text(size=14))+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("Other"='gray70',
                             "Wet Sclerophyll"='darkgreen',
                             "Eucalypt Woodlands" = 'tan3',
                             "Dry Sclerophyll"= 'chartreuse3'))
  
  
p1


#vesta1
p1<-ggplot(df2, aes(x=rate/1000, y=ros_vesta_afdrs/1000, colour=mvg))+
  geom_line(data=loa, aes(x=seq(0, 40, 0.01), y=seq(0, 40, 0.01)), size=0.8, linetype='solid', color='black')+
  geom_line(data=loa, aes(x=seq(0, 90, 0.01)*1.35, y=seq(0, 90, 0.01)), size=0.8, linetype='longdash', color='black')+
  geom_line(data=loa, aes(x=seq(0, 90, 0.01), y=seq(0, 90, 0.01)*1.35), size=0.8, linetype='longdash', color='black')+
  geom_point(size=4)+
  theme_bw()+
  theme(panel.border=element_blank())+
  theme(panel.grid=element_blank())+
  theme(panel.border=element_rect(colour='black', fill='transparent'))+
  theme(panel.grid=element_blank())+
  scale_y_continuous(limits=c(0, 11), breaks = seq(0, 11, by = 2), expand = c(0, 0))+
  scale_x_continuous(limits=c(0, 11), breaks = seq(0, 11, by = 2), expand = c(0, 0))+
  labs(x=expression(Observed~rate~of~spread~(km~hr^-1)), 
       y=expression(Predicted~rate~of~spread~(km~hr^-1)))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=14))+
  theme(axis.text.y= element_text(size=14))+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("Other"='gray70',
                               "Wet Sclerophyll"='darkgreen',
                               "Eucalypt Woodlands" = 'tan3',
                               "Dry Sclerophyll"= 'chartreuse3'))


p1


#10%ws
p1<-ggplot(df2, aes(x=rate/1000, y=ros_10wind/1000, colour=mvg))+
  geom_line(data=loa, aes(x=seq(0, 40, 0.01), y=seq(0, 40, 0.01)), size=0.8, linetype='solid', color='black')+
  geom_line(data=loa, aes(x=seq(0, 90, 0.01)*1.35, y=seq(0, 90, 0.01)), size=0.8, linetype='longdash', color='black')+
  geom_line(data=loa, aes(x=seq(0, 90, 0.01), y=seq(0, 90, 0.01)*1.35), size=0.8, linetype='longdash', color='black')+
  geom_point(size=4)+
  theme_bw()+
  theme(panel.border=element_blank())+
  theme(panel.grid=element_blank())+
  theme(panel.border=element_rect(colour='black', fill='transparent'))+
  theme(panel.grid=element_blank())+
  scale_y_continuous(limits=c(0, 11), breaks = seq(0, 11, by = 2), expand = c(0, 0))+
  scale_x_continuous(limits=c(0, 11), breaks = seq(0, 11, by = 2), expand = c(0, 0))+
  labs(x=expression(Observed~rate~of~spread~(km~hr^-1)), 
       y=expression(Predicted~rate~of~spread~(km~hr^-1)))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=14))+
  theme(axis.text.y= element_text(size=14))+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("Other"='gray70',
                               "Wet Sclerophyll"='darkgreen',
                               "Eucalypt Woodlands" = 'tan3',
                               "Dry Sclerophyll"= 'chartreuse3'))


p1


#mk5
p1<-ggplot(df2, aes(x=rate/1000, y=ros_mk5/1000, colour=mvg))+
  geom_line(data=loa, aes(x=seq(0, 40, 0.01), y=seq(0, 40, 0.01)), size=0.8, linetype='solid', color='black')+
  geom_line(data=loa, aes(x=seq(0, 90, 0.01)*1.35, y=seq(0, 90, 0.01)), size=0.8, linetype='longdash', color='black')+
  geom_line(data=loa, aes(x=seq(0, 90, 0.01), y=seq(0, 90, 0.01)*1.35), size=0.8, linetype='longdash', color='black')+
  geom_point(size=4)+
  theme_bw()+
  theme(panel.border=element_blank())+
  theme(panel.grid=element_blank())+
  theme(panel.border=element_rect(colour='black', fill='transparent'))+
  theme(panel.grid=element_blank())+
  scale_y_continuous(limits=c(0, 11), breaks = seq(0, 11, by = 2), expand = c(0, 0))+
  scale_x_continuous(limits=c(0, 11), breaks = seq(0, 11, by = 2), expand = c(0, 0))+
  labs(x=expression(Observed~rate~of~spread~(km~hr^-1)), 
       y=expression(Predicted~rate~of~spread~(km~hr^-1)))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=14))+
  theme(axis.text.y= element_text(size=14))+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("Other"='gray70',
                               "Wet Sclerophyll"='darkgreen',
                               "Eucalypt Woodlands" = 'tan3',
                               "Dry Sclerophyll"= 'chartreuse3'))


p1

#####

#residual plots

loa2=data.frame()

#vesta2
p1<-ggplot(df2, aes(x=rate/1000, y=res_vesta2/1000, colour=mvg))+
  geom_line(data=loa2, aes(x=seq(0, 40, 0.01), y=seq(0, 0, 0)), size=0.5, linetype='solid', color='black')+
  geom_point(size=4)+
  theme_bw()+
  theme(panel.border=element_blank())+
  theme(panel.grid=element_blank())+
  theme(panel.border=element_rect(colour='black', fill='transparent'))+
  theme(panel.grid=element_blank())+
  scale_y_continuous(limits=c(-6, 5), breaks = seq(-6, 11, by = 2), expand = c(0, 0))+
  scale_x_continuous(limits=c(0, 11), breaks = seq(0, 11, by = 2), expand = c(0, 0))+
  labs(x=expression(Observed~rate~of~spread~(km~hr^-1)), 
       y=expression(Residuals~(km~hr^-1)))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=14))+
  theme(axis.text.y= element_text(size=14))+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("Other"='gray70',
                               "Wet Sclerophyll"='darkgreen',
                               "Dry Sclerophyll"= 'chartreuse3',
                               "Eucalypt Woodlands" = 'tan3'))


p1


#vesta1
p1<-ggplot(df2, aes(x=rate/1000, y=res_vesta1/1000, colour=mvg))+
  geom_line(data=loa2, aes(x=seq(0, 40, 0.01), y=seq(0, 0, 0)), size=0.5, linetype='solid', color='black')+
  geom_point(size=4)+
  theme_bw()+
  theme(panel.border=element_blank())+
  theme(panel.grid=element_blank())+
  theme(panel.border=element_rect(colour='black', fill='transparent'))+
  theme(panel.grid=element_blank())+
  scale_y_continuous(limits=c(-6, 5), breaks = seq(-6, 11, by = 2), expand = c(0, 0))+
  scale_x_continuous(limits=c(0, 11), breaks = seq(0, 11, by = 2), expand = c(0, 0))+
  labs(x=expression(Observed~rate~of~spread~(km~hr^-1)), 
       y=expression(Residuals~(km~hr^-1)))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=14))+
  theme(axis.text.y= element_text(size=14))+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("Other"='gray70',
                               "Wet Sclerophyll"='darkgreen',
                               "Dry Sclerophyll"= 'chartreuse3',
                               "Eucalypt Woodlands" = 'tan3'))


p1


#10wind
p1<-ggplot(df2, aes(x=rate/1000, y=res_10wind/1000, colour=mvg))+
  geom_line(data=loa2, aes(x=seq(0, 40, 0.01), y=seq(0, 0, 0)), size=0.5, linetype='solid', color='black')+
  geom_point(size=4)+
  theme_bw()+
  theme(panel.border=element_blank())+
  theme(panel.grid=element_blank())+
  theme(panel.border=element_rect(colour='black', fill='transparent'))+
  theme(panel.grid=element_blank())+
  scale_y_continuous(limits=c(-6, 5), breaks = seq(-6, 11, by = 2), expand = c(0, 0))+
  scale_x_continuous(limits=c(0, 11), breaks = seq(0, 11, by = 2), expand = c(0, 0))+
  labs(x=expression(Observed~rate~of~spread~(km~hr^-1)), 
       y=expression(Residuals~(km~hr^-1)))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=14))+
  theme(axis.text.y= element_text(size=14))+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("Other"='gray70',
                               "Wet Sclerophyll"='darkgreen',
                               "Dry Sclerophyll"= 'chartreuse3',
                               "Eucalypt Woodlands" = 'tan3'))


p1


#mk5
p1<-ggplot(df2, aes(x=rate/1000, y=res_mk5/1000, colour=mvg))+
  geom_line(data=loa2, aes(x=seq(0, 40, 0.01), y=seq(0, 0, 0)), size=0.5, linetype='solid', color='black')+
  geom_point(size=4)+
  theme_bw()+
  theme(panel.border=element_blank())+
  theme(panel.grid=element_blank())+
  theme(panel.border=element_rect(colour='black', fill='transparent'))+
  theme(panel.grid=element_blank())+
  scale_y_continuous(limits=c(-10, 2), breaks = seq(-6, 11, by = 2), expand = c(0, 0))+
  scale_x_continuous(limits=c(0, 11), breaks = seq(0, 11, by = 2), expand = c(0, 0))+
  labs(x=expression(Observed~rate~of~spread~(km~hr^-1)), 
       y=expression(Residuals~(km~hr^-1)))+
  theme(axis.title.x = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))+
  theme(axis.text.x = element_text(size=14))+
  theme(axis.text.y= element_text(size=14))+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("Other"='gray70',
                               "Wet Sclerophyll"='darkgreen',
                               "Dry Sclerophyll"= 'chartreuse3',
                               "Eucalypt Woodlands" = 'tan3'))


p1


#####

#residual interactions
#RF models with PD interaction plots

rf1 <- randomForest(res_10wind ~ rh + ws_net, data=df2,
                        ntree=500, mtry=2, na.action=na.omit, importance=TRUE)
rf1

pd1 <- partial(rf1, pred.var=c("rh", "ws_net"), progress=TRUE, smooth=TRUE, trim.outliers = TRUE)

plotPartial(pd1, contour=FALSE, at = c(-Inf, seq(-4000, 4000, by = 500), Inf), col.region=grDevices::hcl.colors(100, 'purple-brown', rev=TRUE),
            contour.color='black', lwd=2, xlab=list(label='Relative humidity (%)', fontsize=14),
            ylab=list(label=expression(Windspeed~(km~hr^-1)), fontsize=14),
            smooth=TRUE, smooth.method='gam',
            rug=TRUE, train=df2,
            scales=list(y = list(tck=c(1, 0), at= c(-Inf, seq(0, 80, by = 10), Inf)), 
                        x = list(tck=c(1, 0), at= c(-Inf, seq(0, 100, by = 10), Inf))))


#

rf2 <- randomForest(res_mk5 ~ rh + ws_net, data=df2,
                    ntree=500, mtry=2, na.action=na.omit, importance=TRUE)
rf2

pd2 <- partial(rf2, pred.var=c("rh", "ws_net"), progress=TRUE, smooth=TRUE, trim.outliers = TRUE)

plotPartial(pd2, contour=FALSE, at = c(-Inf, seq(-4000, 4000, by = 500), Inf), col.region=grDevices::hcl.colors(100, 'purple-brown', rev=TRUE),
            contour.color='black', lwd=2, xlab=list(label='Relative humidity (%)', fontsize=14),
            ylab=list(label=expression(Windspeed~(km~hr^-1)), fontsize=14),
            smooth=TRUE, smooth.method='gam',
            rug=TRUE, train=df2,
            scales=list(y = list(tck=c(1, 0), at= c(-Inf, seq(0, 80, by = 10), Inf)), 
                        x = list(tck=c(1, 0), at= c(-Inf, seq(0, 100, by = 10), Inf))))

#

rf3 <- randomForest(res_vesta2 ~ rh + ws_net, data=df2,
                    ntree=500, mtry=2, na.action=na.omit, importance=TRUE)
rf3

pd3 <- partial(rf3, pred.var=c("rh", "ws_net"), progress=TRUE, smooth=TRUE, trim.outliers = TRUE)

plotPartial(pd3, contour=FALSE, at = c(-Inf, seq(-4000, 4000, by = 500), Inf), col.region=grDevices::hcl.colors(100, 'purple-brown', rev=TRUE),
            contour.color='black', lwd=2, xlab=list(label='Relative humidity (%)', fontsize=14),
            ylab=list(label=expression(Windspeed~(km~hr^-1)), fontsize=14),
            rug=TRUE, train=df2,
            scales=list(y = list(tck=c(1, 0), at= c(-Inf, seq(0, 80, by = 10), Inf)), 
                        x = list(tck=c(1, 0), at= c(-Inf, seq(0, 100, by = 10), Inf))))

#



rf4 <- randomForest(res_vesta1 ~ rh + ws_net, data=df2,
                    ntree=500, mtry=2, na.action=na.omit, importance=TRUE)
rf4

pd4 <- partial(rf4, pred.var=c("rh", "ws_net"), progress=TRUE, smooth=TRUE, trim.outliers = TRUE)

plotPartial(pd4, contour=FALSE, at = c(-Inf, seq(-4000, 4000, by = 500), Inf), col.region=grDevices::hcl.colors(100, 'purple-brown', rev=TRUE),
            contour.color='black', lwd=2, xlab=list(label='Relative humidity (%)', fontsize=14),
            ylab=list(label=expression(Windspeed~(km~hr^-1)), fontsize=14),
            smooth=TRUE, smooth.method='gam',
            rug=TRUE, train=df2,
            scales=list(y = list(tck=c(1, 0), at= c(-Inf, seq(0, 80, by = 10), Inf)), 
                        x = list(tck=c(1, 0), at= c(-Inf, seq(0, 100, by = 10), Inf))))



#####




