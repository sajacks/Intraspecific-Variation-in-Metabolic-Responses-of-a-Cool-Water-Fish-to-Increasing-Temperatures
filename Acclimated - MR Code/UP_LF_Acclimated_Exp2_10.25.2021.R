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
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_10.24.2021")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used
info_21.25<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                        Mass=c(15.1,18.51,17.85,16.08,13.67,15.6,14.22,14.16), 
                        Volume = c(317,317,317,317,317,317,317,317), 
                        DO.unit = "mg/L")
info_21.26<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(14.4,18.29,17.49,15.61,13.2,15.38,14.15,13.38), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_23.27<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(18.56,19.14,15.12,20.56,14.96,14.85,17.69,16.52), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_23.28<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(18.4,18.75,14.92,19.98,15.14,14.76,17.68,16.30), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_25.29<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(15.15,17.66,14.82,17.08,18.09,18.23,18.49,14.87), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_25.30<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(13.60,17.41,14.91,16.85,17.93,17.47,18.18,14.62), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")

#Convert blank runs used for background respiration
convert.rMR('empty_pre_10.24.2021_raw.txt', 
            'Converted Empty_Pre 10.24.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_pre_10.26.2021_raw.txt',
            'Converted Empty_Post 10.26.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_pre_10.28.2021_raw.txt',
            'Converted Empty_Post 10.28.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('empty_post_pre_10.30.2021_raw.txt',
            'Converted Empty_post_pre 10.30.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")



#Import Background Respiration Files
LF21.pre<-import.test('Converted Empty_Pre 10.24.2021.txt', 
                      info.data=info_21.3, 
                      logger= "AutoResp",
                      n.chamber=8,
                      plot.temperature= F,
                      plot.oxygen = F)
LF21.post<-import.test('Converted Empty_Post 10.26.2021.txt',
                       info.data=info_21.4, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)
LF23.pre<-import.test('Converted Empty_Post 10.26.2021.txt',
                      info.data=info_23.5, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
LF23.post<-import.test('Converted Empty_Post 10.28.2021.txt',
                       info.data=info_23.6, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)

LF25.pre<-import.test('Converted Empty_Post 10.28.2021.txt',
                      info.data=info_25.7, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
LF25.post<-import.test('Converted Empty_post_pre 10.30.2021.txt',
                       info.data=info_25.8, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)


#######################################################################
#Standard Metabolic Rate
#######################################################################
#Convert %Air Saturation to mg/L for SMR measurements and Empty chamber runs
#First Run of Acute SMR Lower Pennisula from 2019 Fall fingerling data
convert.rMR('SMR_21_LF_10.25.2021_raw.txt',
            'Converted SMR 21 10.25.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_23_LF_10.27.2021_raw.txt',
            'Converted SMR 23 10.27.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_25_LF_10.29.2021_raw.txt',
            'Converted SMR 25 10.29.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw SMR Data
LF.SMR.21.raw<-import.meas('Converted SMR 21 10.25.2021.txt', 
                              info.data=info_21.25, 
                              logger="AutoResp", 
                              n.chamber=8,
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
LF.SMR.23.raw<-import.meas('Converted SMR 23 10.27.2021.txt', 
                              info.data=info_23.27, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
LF.SMR.25.raw<-import.meas('Converted SMR 25 10.29.2021.txt', 
                              info.data=info_25.29, 
                              logger="AutoResp", 
                              n.chamber=8,
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)


#Clean Raw SMR data with background respiration measurments
LF.SMR.21.clean<-correct.meas(info.data=info_21.25,
                              pre.data=LF21.pre,
                              post.data=LF21.post,
                              meas.data=LF.SMR.21.raw,
                              method="exponential")
LF.SMR.23.clean<-correct.meas(info.data=info_23.27,
                              pre.data=LF23.pre,
                              post.data=LF23.post, 
                              meas.data=LF.SMR.23.raw,
                              method="exponential")
LF.SMR.25.clean<-correct.meas(info.data=info_25.29,
                              pre.data=LF25.pre,
                              post.data=LF25.post,
                              meas.data=LF.SMR.25.raw,
                              method="exponential")

#Extract slopes for Standard Metabolic Rate
LF.SMR.21.slopes<-extract.slope(LF.SMR.21.clean,
                                 method="calcSMR.quant",
                                 r2=0.9,
                                 p=0.2)
LF.SMR.23.slopes<-extract.slope(LF.SMR.23.clean,
                                 method="calcSMR.quant",
                                 r2=c(0.9),
                                 p=0.2)
LF.SMR.25.slopes<-extract.slope(LF.SMR.25.clean,
                                 method="calcSMR.quant",
                                 r2=c(0.9),
                                 p=0.2)


##Extract slopes for Maximum Metabolic Rate during SMR Measurements
LF.SMR.MMR.21.slopes<-extract.slope(LF.SMR.21.clean,
                                     method="max",
                                     n.slope=1)
LF.SMR.MMR.23.slopes<-extract.slope(LF.SMR.23.clean,
                                     method="max",
                                     n.slope=1)
LF.SMR.MMR.25.slopes<-extract.slope(LF.SMR.25.clean,
                                       method="max",
                                       n.slope=1)



##Calculate Standard Metabolic Rates
LF.SMR.21<-calculate.MR(LF.SMR.21.slopes,density=1000)
LF.SMR.23<-calculate.MR(LF.SMR.23.slopes,density=1000)
LF.SMR.25<-calculate.MR(LF.SMR.25.slopes,density=1000)


##Calculate Maximum Metabolic Rates 
LF.SMR.MMR.21<-calculate.MR(LF.SMR.MMR.21.slopes,density=1000)
LF.SMR.MMR.23<-calculate.MR(LF.SMR.MMR.23.slopes,density=1000)
LF.SMR.MMR.25<-calculate.MR(LF.SMR.MMR.25.slopes,density=1000)


##Combine Data into one Dataframe
LF.SMR.Acclimated<-rbind.data.frame(LF.SMR.21,LF.SMR.23,LF.SMR.25)
LF.SMR.MMR.Acclimated<-rbind.data.frame(LF.SMR.MMR.21,LF.SMR.MMR.23,LF.SMR.MMR.25)

##Add additional Columns
Temperature<-c(21,21,21,21,21,21,21,21,23,23,23,23,23,23,23,23,25,25,25,25,25,25,25,25)
Trial<-c("SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR")
Event<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3)

LF.SMR.Acclimated<-cbind(LF.SMR.Acclimated,Trial,Temperature,Event)
LF.SMR.MMR.Acclimated<-cbind(LF.SMR.MMR.Acclimated,Trial,Temperature,Event)

##############33
##Save out CSVs


write.csv(LF.SMR.Acclimated,"LF Acclimated SMR EXP 2.csv")
write.csv(LF.SMR.MMR.Acclimated,"LF Acclimated SMR-MMR EXP 2.csv")


#######################################################################
#Maximum Metabolic Rate
########################################################################

#Convert %Air Saturation to mg/L for MMR measurements and Empty chamber runs
#First Run of Acute MMR Upper Pennisula
convert.rMR('MMR_21_LF_10.26.2021_raw.txt',
            'Converted MMR 21 10.26.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_23_LF_10.28.2021_raw.txt',
            'Converted MMR 23 10.28.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_25_LF_10.30.2021_raw.txt',
            'Converted MMR 25 10.30.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw MMR Data
LF.MMR.21.raw<-import.meas('Converted MMR 21 10.26.2021.txt', 
                            info.data=info_21.26, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
LF.MMR.23.raw<-import.meas('Converted MMR 23 10.28.2021.txt', 
                            info.data=info_23.28, 
                            n.chamber=8,
                            logger="AutoResp", 
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
LF.MMR.25.raw<-import.meas('Converted MMR 25 10.30.2021.txt', 
                            info.data=info_25.30, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)

#Clean Raw MMR data with background respiration measurments
LF.MMR.21.clean<-correct.meas(info.data=info_21.26,
                              pre.data=LF21.pre,
                              post.data=LF21.post,
                              meas.data=LF.MMR.21.raw,
                              method="exponential")
LF.MMR.23.clean<-correct.meas(info.data=info_23.28,
                              pre.data=LF23.pre,
                              post.data=LF23.post, 
                              meas.data=LF.MMR.23.raw,
                              method="exponential")
LF.MMR.25.clean<-correct.meas(info.data=info_25.30,
                              pre.data=LF25.pre,
                              post.data=LF25.post,
                              meas.data=LF.MMR.25.raw,
                              method="exponential")

#Extract slopes
LF.MMR.21.slopes<-extract.slope(LF.MMR.21.clean,
                                 method="max",
                                 r2=0.95,
                                 n.slope=1)
LF.MMR.23.slopes<-extract.slope(LF.MMR.23.clean,
                                 method="max",
                                 r2=0.95,
                                 n.slope=1)
LF.MMR.25.slopes<-extract.slope(LF.MMR.25.clean,
                                 method="max",
                                 r2=0.95,
                                 n.slope=1)

##Calculate Metabolic Rates
LF.MMR.21<-calculate.MR(LF.MMR.21.slopes,density=1000)
LF.MMR.23<-calculate.MR(LF.MMR.23.slopes,density=1000)
LF.MMR.25<-calculate.MR(LF.MMR.25.slopes,density=1000)


##Combine data into one Dataframe
LF.MMR.Acclimated<-rbind.data.frame(LF.MMR.21,LF.MMR.23,LF.MMR.25)

##Add additional columns
Trial<-c("MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR")
Temperature<-c(21,21,21,21,21,21,21,21,23,23,23,23,23,23,23,23,25,25,25,25,25,25,25,25)
Event<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3)
LF.MMR.Acclimated<-cbind(LF.MMR.Acclimated,Trial,Temperature,Event)

##Save data our as CSV


write.csv(LF.MMR.Acclimated,"LP Acclimated MMR EXP 2.csv")

#######################################################################
#Data Analysis and grapahs
#######################################################################

LF.MMR<-rbind(LF.SMR.MMR.Acclimated,LF.MMR.Acclimated)

grid.arrange(
  ggplot(subset(LF.MMR,Event%in%c(1)), 
         aes(x = Trial, y = MR.mass, group=Ind,color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, UP.1 21") + theme(legend.position="none")+ylim(200,600),
  ggplot(subset(LF.MMR,Event%in%c(2)), 
         aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, UP.1 23") + theme(legend.position="none")+ylim(200,600),
  ggplot(subset(LF.MMR,Event%in%c(3)), 
         aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, UP.1 25") + theme()+ylim(200,600),
  
  ncol=3
)
