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
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_10.10.2021")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used
info_21.10<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                        Mass=c(10.72,9.97,11.57,13.67,13.26,16.34,16.76,10.62), 
                        Volume = c(317,317,317,317,317,317,317,317), 
                        DO.unit = "mg/L")
info_21.11<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(10.65,9.78,11.10,13.55,13.18,16.25,16.63,10.64), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_23.12<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(11.91,17.75,11.09,12.65,9.38,17.61,8.79,12.49), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_23.13<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(11.69,17.15,10.99,12.12,9.2,17.30,8.58,12.28), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_25.14<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(14.08,16.23,11.92,16.58, 14.47,11.97,11.63,17.31), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_25.15<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(13.86,15.92,11.55,16.18,14.36,11.63,11.28,16.85), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")

#Convert blank runs used for background respiration
convert.rMR('empty_pre_10.9.2021_raw.txt', 
            'Converted Empty_Pre 10.09.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_pre_10.11.2021_raw.txt',
            'Converted Empty_Post 10.11.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_pre_10.13.2021_raw.txt',
            'Converted Empty_Post 10.13.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('empty_post_10.15.2021_raw.txt',
            'Converted Empty_post 10.15.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")



#Import Background Respiration Files
BC21.pre<-import.test('Converted Empty_Pre 10.09.2021.txt', 
                      info.data=info_21.10, 
                      logger= "AutoResp",
                      n.chamber=8,
                      plot.temperature= F,
                      plot.oxygen = F)
BC21.post<-import.test('Converted Empty_Post 10.11.2021.txt',
                       info.data=info_21.11, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)
BC23.pre<-import.test('Converted Empty_Post 10.11.2021.txt',
                      info.data=info_23.12, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
BC23.post<-import.test('Converted Empty_Post 10.13.2021.txt',
                       info.data=info_23.13, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)

BC25.pre<-import.test('Converted Empty_Post 10.13.2021.txt',
                      info.data=info_25.14, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
BC25.post<-import.test('Converted Empty_post 10.15.2021.txt',
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
convert.rMR('SMR_21_BC_10.10.2021_raw.txt',
            'Converted SMR 21 10.10.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_23_BC_10.12.2021_raw.txt',
            'Converted SMR 23 10.12.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_25_BC_10.14.2021_raw.txt',
            'Converted SMR 25 10.14.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw SMR Data
BC.SMR.21.raw<-import.meas('Converted SMR 21 10.10.2021.txt', 
                              info.data=info_21.10, 
                              logger="AutoResp", 
                              n.chamber=8,
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
BC.SMR.23.raw<-import.meas('Converted SMR 23 10.12.2021.txt', 
                              info.data=info_23.12, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
BC.SMR.25.raw<-import.meas('Converted SMR 25 10.14.2021.txt', 
                              info.data=info_25.14, 
                              logger="AutoResp", 
                              n.chamber=8,
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)


#Clean Raw SMR data with background respiration measurments
BC.SMR.21.clean<-correct.meas(info.data=info_21.10,
                              pre.data=BC21.pre,
                              post.data=BC21.post,
                              meas.data=BC.SMR.21.raw,
                              method="exponential")
BC.SMR.23.clean<-correct.meas(info.data=info_23.12,
                              pre.data=BC23.pre,
                              post.data=BC23.post, 
                              meas.data=BC.SMR.23.raw,
                              method="exponential")
BC.SMR.25.clean<-correct.meas(info.data=info_25.14,
                              pre.data=BC25.pre,
                              post.data=BC25.post,
                              meas.data=BC.SMR.25.raw,
                              method="exponential")

#Extract slopes for Standard Metabolic Rate
BC.SMR.21.slopes<-extract.slope(BC.SMR.21.clean,
                                 method="calcSMR.quant",
                                 r2=0.9,
                                 p=0.2)
BC.SMR.23.slopes<-extract.slope(BC.SMR.23.clean,
                                 method="calcSMR.quant",
                                 r2=c(0.9),
                                 p=0.2)
BC.SMR.25.slopes<-extract.slope(BC.SMR.25.clean,
                                 method="calcSMR.quant",
                                 r2=c(0.9),
                                 p=0.2)


##Extract slopes for Maximum Metabolic Rate during SMR Measurements
BC.SMR.MMR.21.slopes<-extract.slope(BC.SMR.21.clean,
                                     method="max",
                                     n.slope=1)
BC.SMR.MMR.23.slopes<-extract.slope(BC.SMR.23.clean,
                                     method="max",
                                     n.slope=1)
BC.SMR.MMR.25.slopes<-extract.slope(BC.SMR.25.clean,
                                       method="max",
                                       n.slope=1)



##Calculate Standard Metabolic Rates
BC.SMR.21<-calculate.MR(BC.SMR.21.slopes,density=1000)
BC.SMR.23<-calculate.MR(BC.SMR.23.slopes,density=1000)
BC.SMR.25<-calculate.MR(BC.SMR.25.slopes,density=1000)


##Calculate Maximum Metabolic Rates 
BC.SMR.MMR.21<-calculate.MR(BC.SMR.MMR.21.slopes,density=1000)
BC.SMR.MMR.23<-calculate.MR(BC.SMR.MMR.23.slopes,density=1000)
BC.SMR.MMR.25<-calculate.MR(BC.SMR.MMR.25.slopes,density=1000)


##Combine Data into one Dataframe
BC.SMR.Acclimated<-rbind.data.frame(BC.SMR.21,BC.SMR.23,BC.SMR.25)
BC.SMR.MMR.Acclimated<-rbind.data.frame(BC.SMR.MMR.21,BC.SMR.MMR.23,BC.SMR.MMR.25)

##Add additional Columns
Temperature<-c(21,21,21,21,21,21,21,21,23,23,23,23,23,23,23,23,25,25,25,25,25,25,25,25)
Trial<-c("SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR")
Event<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3)

BC.SMR.Acclimated<-cbind(BC.SMR.Acclimated,Trial,Temperature,Event)
BC.SMR.MMR.Acclimated<-cbind(BC.SMR.MMR.Acclimated,Trial,Temperature,Event)

##############33
##Save out CSVs


write.csv(BC.SMR.Acclimated,"BC Acclimated SMR EXP 1.csv")
write.csv(BC.SMR.MMR.Acclimated,"BC Acclimated SMR-MMR EXP 1.csv")


#######################################################################
#Maximum Metabolic Rate
########################################################################

#Convert %Air Saturation to mg/L for MMR measurements and Empty chamber runs
#First Run of Acute MMR Upper Pennisula
convert.rMR('MMR_21_BC_10.11.2021_raw.txt',
            'Converted MMR 21 10.11.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_23_BC_10.13.2021_raw.txt',
            'Converted MMR 23 10.13.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_25_BC_10.15.2021_raw.txt',
            'Converted MMR 25 10.15.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw MMR Data
BC.MMR.21.raw<-import.meas('Converted MMR 21 10.11.2021.txt', 
                            info.data=info_21.11, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
BC.MMR.23.raw<-import.meas('Converted MMR 23 10.13.2021.txt', 
                            info.data=info_23.13, 
                            n.chamber=8,
                            logger="AutoResp", 
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
BC.MMR.25.raw<-import.meas('Converted MMR 25 10.15.2021.txt', 
                            info.data=info_25.15, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)

#Clean Raw MMR data with background respiration measurments
BC.MMR.21.clean<-correct.meas(info.data=info_21.11,
                              pre.data=BC21.pre,
                              post.data=BC21.post,
                              meas.data=BC.MMR.21.raw,
                              method="exponential")
BC.MMR.23.clean<-correct.meas(info.data=info_23.13,
                              pre.data=BC23.pre,
                              post.data=BC23.post, 
                              meas.data=BC.MMR.23.raw,
                              method="exponential")
BC.MMR.25.clean<-correct.meas(info.data=info_25.15,
                              pre.data=BC25.pre,
                              post.data=BC25.post,
                              meas.data=BC.MMR.25.raw,
                              method="exponential")

#Extract slopes
BC.MMR.21.slopes<-extract.slope(BC.MMR.21.clean,
                                 method="max",
                                 r2=0.95,
                                 n.slope=1)
BC.MMR.23.slopes<-extract.slope(BC.MMR.23.clean,
                                 method="max",
                                 r2=0.95,
                                 n.slope=1)
BC.MMR.25.slopes<-extract.slope(BC.MMR.25.clean,
                                 method="max",
                                 r2=0.95,
                                 n.slope=1)

##Calculate Metabolic Rates
BC.MMR.21<-calculate.MR(BC.MMR.21.slopes,density=1000)
BC.MMR.23<-calculate.MR(BC.MMR.23.slopes,density=1000)
BC.MMR.25<-calculate.MR(BC.MMR.25.slopes,density=1000)


##Combine data into one Dataframe
BC.MMR.Acclimated<-rbind.data.frame(BC.MMR.21,BC.MMR.23,BC.MMR.25)

##Add additional columns
Trial<-c("MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR")
Temperature<-c(21,21,21,21,21,21,21,21,23,23,23,23,23,23,23,23,25,25,25,25,25,25,25,25)
Event<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3)
BC.MMR.Acclimated<-cbind(BC.MMR.Acclimated,Trial,Temperature,Event)

##Save data our as CSV


write.csv(BC.MMR.Acclimated,"BC Acclimated MMR EXP 1.csv")

#######################################################################
#Data Analysis and grapahs
#######################################################################

BC.MMR<-rbind(BC.SMR.MMR.Acclimated,BC.MMR.Acclimated)

grid.arrange(
  ggplot(subset(BC.MMR,Event%in%c(1)), 
         aes(x = Trial, y = MR.mass, group=Ind,color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, UP.1 21") + theme(legend.position="none")+ylim(200,600),
  ggplot(subset(BC.MMR,Event%in%c(2)), 
         aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, UP.1 23") + theme(legend.position="none")+ylim(200,600),
  ggplot(subset(BC.MMR,Event%in%c(3)), 
         aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, UP.1 25") + theme(legend.position="none")+ylim(200,600),
  
  ncol=3
)
