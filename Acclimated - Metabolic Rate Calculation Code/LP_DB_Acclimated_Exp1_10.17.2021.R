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
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_10.17.2021")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used
info_21.17<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                        Mass=c(13.87,9.23,12.37,12.34,10.220,14.066,13.914,13.9975), 
                        Volume = c(317,317,317,317,317,317,317,317), 
                        DO.unit = "mg/L")
info_21.18<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(13.82,8.92,12.11,12.21,10.05,13.64,13.77,13.98), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_23.19<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(10.51,13.57,10.92,14.26,7.9,12.53, 10.23, 9.91),
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_23.20<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(10.2,13.29,10.67,13.84,7.67,11.94,9.74,9.38), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_25.21<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(13.9,11.99,10.86,12.53,12.09,11.34,14.13,15.67), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_25.22<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(13.32,11.79,10.55,12.18,11.84,11.26,14.13,15.67), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")

#Convert blank runs used for background respiration
convert.rMR('empty_pre_10.16.2021_raw.txt', 
            'Converted Empty_Pre 10.16.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_pre_10.18.2021_raw.txt',
            'Converted Empty_Post 10.18.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_pre_10.20.2021_raw.txt',
            'Converted Empty_Post 10.20.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('empty_post_10.22.2021_raw.txt',
            'Converted Empty_post 10.22.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")



#Import Background Respiration Files
DB21.pre<-import.test('Converted Empty_Pre 10.16.2021.txt', 
                      info.data=info_21.10, 
                      logger= "AutoResp",
                      n.chamber=8,
                      plot.temperature= F,
                      plot.oxygen = F)
DB21.post<-import.test('Converted Empty_Post 10.18.2021.txt',
                       info.data=info_21.11, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)
DB23.pre<-import.test('Converted Empty_Post 10.18.2021.txt',
                      info.data=info_23.12, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
DB23.post<-import.test('Converted Empty_Post 10.20.2021.txt',
                       info.data=info_23.13, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)

DB25.pre<-import.test('Converted Empty_Post 10.20.2021.txt',
                      info.data=info_25.14, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
DB25.post<-import.test('Converted Empty_post 10.22.2021.txt',
                       info.data=info_25.15, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)


#######################################################################
#Standard Metabolic Rate
#######################################################################
#Convert %Air Saturation to mg/L for SMR measurements and Empty chamber runs
#First Run of Acute SMR Lower Pennisula from 2019 Fall fingerling data
convert.rMR('SMR_21_DB_10.17.2021_raw.txt',
            'Converted SMR 21 10.17.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_23_DB_10.19.2021_raw.txt',
            'Converted SMR 23 10.19.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_25_DB_10.21.2021_raw.txt',
            'Converted SMR 25 10.21.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw SMR Data
DB.SMR.21.raw<-import.meas('Converted SMR 21 10.17.2021.txt', 
                              info.data=info_21.17, 
                              logger="AutoResp", 
                              n.chamber=8,
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
DB.SMR.23.raw<-import.meas('Converted SMR 23 10.19.2021.txt', 
                              info.data=info_23.19, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
DB.SMR.25.raw<-import.meas('Converted SMR 25 10.21.2021.txt', 
                              info.data=info_25.21, 
                              logger="AutoResp", 
                              n.chamber=8,
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)


#Clean Raw SMR data with background respiration measurments
DB.SMR.21.clean<-correct.meas(info.data=info_21.17,
                              pre.data=DB21.pre,
                              post.data=DB21.post,
                              meas.data=DB.SMR.21.raw,
                              method="exponential")
DB.SMR.23.clean<-correct.meas(info.data=info_23.19,
                              pre.data=DB23.pre,
                              post.data=DB23.post, 
                              meas.data=DB.SMR.23.raw,
                              method="exponential")
DB.SMR.25.clean<-correct.meas(info.data=info_25.21,
                              pre.data=DB25.pre,
                              post.data=DB25.post,
                              meas.data=DB.SMR.25.raw,
                              method="exponential")

#Extract slopes for Standard Metabolic Rate
DB.SMR.21.slopes<-extract.slope(DB.SMR.21.clean,
                                 method="calcSMR.quant",
                                 r2=0.9,
                                 p=0.2)
DB.SMR.23.slopes<-extract.slope(DB.SMR.23.clean,
                                 method="calcSMR.quant",
                                 r2=c(0.9),
                                 p=0.2)
DB.SMR.25.slopes<-extract.slope(DB.SMR.25.clean,
                                 method="calcSMR.quant",
                                 r2=c(0.9),
                                 p=0.2)


##Extract slopes for Maximum Metabolic Rate during SMR Measurements
DB.SMR.MMR.21.slopes<-extract.slope(DB.SMR.21.clean,
                                     method="max",
                                     n.slope=1)
DB.SMR.MMR.23.slopes<-extract.slope(DB.SMR.23.clean,
                                     method="max",
                                     n.slope=1)
DB.SMR.MMR.25.slopes<-extract.slope(DB.SMR.25.clean,
                                       method="max",
                                       n.slope=1)



##Calculate Standard Metabolic Rates
DB.SMR.21<-calculate.MR(DB.SMR.21.slopes,density=1000)
DB.SMR.23<-calculate.MR(DB.SMR.23.slopes,density=1000)
DB.SMR.25<-calculate.MR(DB.SMR.25.slopes,density=1000)


##Calculate Maximum Metabolic Rates 
DB.SMR.MMR.21<-calculate.MR(DB.SMR.MMR.21.slopes,density=1000)
DB.SMR.MMR.23<-calculate.MR(DB.SMR.MMR.23.slopes,density=1000)
DB.SMR.MMR.25<-calculate.MR(DB.SMR.MMR.25.slopes,density=1000)


##Combine Data into one Dataframe
DB.SMR.Acclimated<-rbind.data.frame(DB.SMR.21,DB.SMR.23,DB.SMR.25)
DB.SMR.MMR.Acclimated<-rbind.data.frame(DB.SMR.MMR.21,DB.SMR.MMR.23,DB.SMR.MMR.25)

##Add additional Columns
Temperature<-c(21,21,21,21,21,21,21,21,23,23,23,23,23,23,23,23,25,25,25,25,25,25,25,25)
Trial<-c("SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR")
Event<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3)

DB.SMR.Acclimated<-cbind(DB.SMR.Acclimated,Trial,Temperature,Event)
DB.SMR.MMR.Acclimated<-cbind(DB.SMR.MMR.Acclimated,Trial,Temperature,Event)

##############33
##Save out CSVs


write.csv(DB.SMR.Acclimated,"DB Acclimated SMR EXP 1.csv")
write.csv(DB.SMR.MMR.Acclimated,"DB Acclimated SMR-MMR EXP 1.csv")


#######################################################################
#Maximum Metabolic Rate
########################################################################

#Convert %Air Saturation to mg/L for MMR measurements and Empty chamber runs
#First Run of Acute MMR Upper Pennisula
convert.rMR('MMR_21_DB_10.18.2021_raw.txt',
            'Converted MMR 21 10.18.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_23_DB_10.20.2021_raw.txt',
            'Converted MMR 23 10.20.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_25_DB_10.22.2021_raw.txt',
            'Converted MMR 25 10.22.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw MMR Data
DB.MMR.21.raw<-import.meas('Converted MMR 21 10.18.2021.txt', 
                            info.data=info_21.18, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
DB.MMR.23.raw<-import.meas('Converted MMR 23 10.20.2021.txt', 
                            info.data=info_23.20, 
                            n.chamber=8,
                            logger="AutoResp", 
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
DB.MMR.25.raw<-import.meas('Converted MMR 25 10.22.2021.txt', 
                            info.data=info_25.22, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)

#Clean Raw MMR data with background respiration measurments
DB.MMR.21.clean<-correct.meas(info.data=info_21.18,
                              pre.data=DB21.pre,
                              post.data=DB21.post,
                              meas.data=DB.MMR.21.raw,
                              method="exponential")
DB.MMR.23.clean<-correct.meas(info.data=info_23.20,
                              pre.data=DB23.pre,
                              post.data=DB23.post, 
                              meas.data=DB.MMR.23.raw,
                              method="exponential")
DB.MMR.25.clean<-correct.meas(info.data=info_25.22,
                              pre.data=DB25.pre,
                              post.data=DB25.post,
                              meas.data=DB.MMR.25.raw,
                              method="exponential")

#Extract slopes
DB.MMR.21.slopes<-extract.slope(DB.MMR.21.clean,
                                 method="max",
                                 r2=0.95,
                                 n.slope=1)
DB.MMR.23.slopes<-extract.slope(DB.MMR.23.clean,
                                 method="max",
                                 r2=0.95,
                                 n.slope=1)
DB.MMR.25.slopes<-extract.slope(DB.MMR.25.clean,
                                 method="max",
                                 r2=0.95,
                                 n.slope=1)

##Calculate Metabolic Rates
DB.MMR.21<-calculate.MR(DB.MMR.21.slopes,density=1000)
DB.MMR.23<-calculate.MR(DB.MMR.23.slopes,density=1000)
DB.MMR.25<-calculate.MR(DB.MMR.25.slopes,density=1000)


##Combine data into one Dataframe
DB.MMR.Acclimated<-rbind.data.frame(DB.MMR.21,DB.MMR.23,DB.MMR.25)

##Add additional columns
Trial<-c("MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR")
Temperature<-c(21,21,21,21,21,21,21,21,23,23,23,23,23,23,23,23,25,25,25,25,25,25,25,25)
Event<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3)
DB.MMR.Acclimated<-cbind(DB.MMR.Acclimated,Trial,Temperature,Event)

##Save data our as CSV


write.csv(DB.MMR.Acclimated,"DB Acclimated MMR EXP 1.csv")

#######################################################################
#Data Analysis and grapahs
#######################################################################

DB.MMR<-rbind(DB.SMR.MMR.Acclimated,DB.MMR.Acclimated)

grid.arrange(
  ggplot(subset(DB.MMR,Event%in%c(1)), 
         aes(x = Trial, y = MR.mass, group=Ind,color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, DB 21") + theme(legend.position="none")+ylim(300,600),
  ggplot(subset(DB.MMR,Event%in%c(2)), 
         aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, DB 23") + theme(legend.position="none")+ylim(300,600),
  ggplot(subset(DB.MMR,Event%in%c(3)), 
         aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, DB 25") + theme(legend.position="none")+ylim(300,600),
  
  ncol=3
)
