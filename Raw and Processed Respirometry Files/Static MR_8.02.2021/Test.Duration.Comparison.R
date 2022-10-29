library(FishResp)
library(tidyverse)
library(knitr)
library(nlme)
library(emmeans)
library(kableExtra)
library(GGally)
library(qqplotr)
library(gridExtra)
library(car)


###Start with BC Experiments 21C
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_8.2.2021")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used

info_21.03<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(4.05,2.89,3.38,3.03,2.71,2.77,2.58,2.72), 
                       Volume = c(78,78,78,78,78,78,78,78), 
                       DO.unit = "mg/L")
info_21.04<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(4.15,2.84,3.13,2.9,2.54,2.71,2.51,2.56), 
                       Volume = c(78,78,78,78,78,78,78,78), 
                       DO.unit = "mg/L")
info_23.05<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(3.65,2.67,3.06,2.75,2.44,2.55,2.53,2.52), 
                       Volume = c(78,78,78,78,78,78,78,78), 
                       DO.unit = "mg/L")
info_25.06<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(3.6,2.69,3.03,2.73,2.37,2.43,2.22,2.51), 
                       Volume = c(78,78,78,78,78,78,78,78), 
                       DO.unit = "mg/L")
info_21.2.07<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                         Mass=c(3.55,2.61,2.99,2.65,2.35,2.44,2.22,2.46), 
                         Volume = c(78,78,78,78,78,78,78,78), 
                         DO.unit = "mg/L")

#Import Background Respiration Files
BC1.pre<-import.test('Converted Empty_Pre 08.03.2021.txt', 
                     info.data=info_21.04, 
                     logger= "AutoResp",
                     n.chamber=8,
                     plot.temperature= F,
                     plot.oxygen = F)
BC1.post<-import.test('Converted Empty_Post 08.07.2021.txt',
                      info.data=info_21.04, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)

#Import Raw SMR Data
BC1.SMR.21.raw<-import.meas('Converted SMR 21 08.03.2021.txt', 
                            info.data=info_21.04, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
BC1.SMR.21.2.raw<-import.meas('Converted SMR 21.2 08.06.2021.txt', 
                              info.data=info_21.2.07, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
BC1.SMR.23.raw<-import.meas('Converted SMR 23 08.04.2021.txt', 
                            info.data=info_23.05, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
BC1.SMR.25.raw<-import.meas('Converted SMR 25 08.05.2021.txt', 
                            info.data=info_25.06, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)

BC1.SMR.21.raw2<-import.meas('Converted SMR 21 08.03.2021.txt', 
                            info.data=info_21.04, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F,
                            start.measure="21:00:00",
                            stop.measure="09:00:00")
BC1.SMR.21.2.raw2<-import.meas('Converted SMR 21.2 08.06.2021.txt', 
                              info.data=info_21.2.07, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F,
                              start.measure="23:27:00",
                              stop.measure="11:27:00")
BC1.SMR.23.raw2<-import.meas('Converted SMR 23 08.04.2021.txt', 
                            info.data=info_23.05, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F,
                            start.measure="19:05:00",
                            stop.measure="07:05:00")
BC1.SMR.25.raw2<-import.meas('Converted SMR 25 08.05.2021.txt', 
                            info.data=info_25.06, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F,
                            start.measure="17:10:00",
                            stop.measure="07:10:00")

#Clean Raw SMR data with background respiration measurments
BC1.SMR.21.clean<-correct.meas(info.data=info_21.04,
                               pre.data=BC1.pre,
                               post.data=BC1.post,
                               meas.data=BC1.SMR.21.raw,
                               method="exponential")
BC1.SMR.21.2.clean<-correct.meas(info.data=info_21.2.07,
                                 pre.data=BC1.pre,
                                 post.data=BC1.post, 
                                 meas.data=BC1.SMR.21.2.raw,
                                 method="exponential")
BC1.SMR.23.clean<-correct.meas(info.data=info_23.05,
                               pre.data=BC1.pre,
                               post.data=BC1.post, 
                               meas.data=BC1.SMR.23.raw,
                               method="exponential")
BC1.SMR.25.clean<-correct.meas(info.data=info_25.06,
                               pre.data=BC1.pre,
                               post.data=BC1.post,
                               meas.data=BC1.SMR.25.raw,
                               method="exponential")

BC1.SMR.21.clean2<-correct.meas(info.data=info_21.04,
                               pre.data=BC1.pre,
                               post.data=BC1.post,
                               meas.data=BC1.SMR.21.raw2,
                               method="exponential")
BC1.SMR.21.2.clean2<-correct.meas(info.data=info_21.2.07,
                                 pre.data=BC1.pre,
                                 post.data=BC1.post, 
                                 meas.data=BC1.SMR.21.2.raw2,
                                 method="exponential")
BC1.SMR.23.clean2<-correct.meas(info.data=info_23.05,
                               pre.data=BC1.pre,
                               post.data=BC1.post, 
                               meas.data=BC1.SMR.23.raw2,
                               method="exponential")
BC1.SMR.25.clean2<-correct.meas(info.data=info_25.06,
                               pre.data=BC1.pre,
                               post.data=BC1.post,
                               meas.data=BC1.SMR.25.raw2,
                               method="exponential")

##Extract slopes for Standard Metabolic Rate
BC1.SMR.21.slopes<-extract.slope(BC1.SMR.21.clean,
                                 method="calcSMR.quant",
                                 r2=0.9,
                                 p=0.2)
BC1.SMR.21.2.slopes<-extract.slope(BC1.SMR.21.2.clean,
                                   method="calcSMR.quant",
                                   r2=0.9,
                                   p=0.2)
BC1.SMR.23.slopes<-extract.slope(BC1.SMR.23.clean,
                                 method="calcSMR.quant",
                                 r2=0.9,
                                 p=0.2)
BC1.SMR.25.slopes<-extract.slope(BC1.SMR.25.clean,
                                 method="calcSMR.quant",
                                 r2=0.9,
                                 p=0.2)

BC1.SMR.21.slopes2<-extract.slope(BC1.SMR.21.clean2,
                                 method="calcSMR.quant",
                                 r2=0.9,
                                 p=0.2)
BC1.SMR.21.2.slopes2<-extract.slope(BC1.SMR.21.2.clean2,
                                   method="calcSMR.quant",
                                   r2=0.9,
                                   p=0.2)
BC1.SMR.23.slopes2<-extract.slope(BC1.SMR.23.clean2,
                                 method="calcSMR.quant",
                                 r2=0.9,
                                 p=0.2)
BC1.SMR.25.slopes2<-extract.slope(BC1.SMR.25.clean2,
                                 method="calcSMR.quant",
                                 r2=0.9,
                                 p=0.2)



##Calculate Standard Metabolic Rates
BC1.SMR.21<-calculate.MR(BC1.SMR.21.slopes,density=1000)
BC1.SMR.21.2<-calculate.MR(BC1.SMR.21.2.slopes,density=1000)
BC1.SMR.23<-calculate.MR(BC1.SMR.23.slopes,density=1000)
BC1.SMR.25<-calculate.MR(BC1.SMR.25.slopes,density=1000)

BC1.SMR.21.12<-calculate.MR(BC1.SMR.21.slopes2,density=1000)
BC1.SMR.21.2.12<-calculate.MR(BC1.SMR.21.2.slopes2,density=1000)
BC1.SMR.23.12<-calculate.MR(BC1.SMR.23.slopes2,density=1000)
BC1.SMR.25.12<-calculate.MR(BC1.SMR.25.slopes2,density=1000)


##Combine Data into one dataframe
BC1.SMR.Acute<-rbind.data.frame(BC1.SMR.21,BC1.SMR.21.2,BC1.SMR.23,BC1.SMR.25)
BC1.SMR.Acute.12h<-rbind.data.frame(BC1.SMR.21.12,BC1.SMR.21.2.12,BC1.SMR.23.12,BC1.SMR.25.12)

##Add additional columns
#Temperature as a Factor
Temperature<-c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,23,23,23,23,23,23,23,23,25,25,25,25,25,25,25,25)
#Trial measurement for SMR.MMR
Trial<-c("SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR","SMR","SMR")
#Experiment Order
Event<-c(1,1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3)

#Trial Length
BC1.SMR.Acute$Duration<-"16"
BC1.SMR.Acute.12h$Duration<-"12"

BC1.SMR.MMR.Acute<-cbind(BC1.SMR.MMR.Acute,Trial,Temperature,Event)
BC1.SMR.Acute<-cbind(BC1.SMR.Acute,Temperature, Event)

BC1.SMR<-rbind(BC1.SMR.Acute,BC1.SMR.Acute.12h)


##########Analysis

ggplot(BC1.SMR,aes(x=Duration,y=MR.mass))+geom_boxplot()+geom_point()








