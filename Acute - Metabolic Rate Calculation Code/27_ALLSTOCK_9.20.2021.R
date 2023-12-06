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
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_9.20.2021")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used
info_27.20<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(17.01,14.51,16.17,12.52,11.74,14.24,11.70,10.05), 
                    Volume = c(317,317,317,317,317,317,317,317), 
                    DO.unit = "mg/L")
info_27.21<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(16.56,14.01,15.39,12.05,11.28,13.48,11.49,9.40), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_27.22<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(8.74,8.75,8.13,8.08,7.93,7.09,10.17,7.98), 
                    Volume = c(317,317,317,317,317,317,317,317), 
                    DO.unit = "mg/L")
info_27.23<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                                Mass=c(8.55,8.46,7.95,8.33,7.54,6.76,9.96,7.72), 
                                Volume = c(317,317,317,317,317,317,317,317), 
                                DO.unit = "mg/L")
info_27.29<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                        Mass=c(9.54,14.82,11.56,12.17,17.98,10.12,13.12,11.69), 
                        Volume = c(317,317,317,317,317,317,317,317), 
                        DO.unit = "mg/L")
info_27.30<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(9.11,14.49,11.07,11.51,17.40,9.33,12.52,11.18), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")

info_31.01<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(19.42,15.68,0,0,12.64,13.05,0,15.20), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_31.02<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(20.43,15.78,0,0,12.73,12.68,0,14.74), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")

info_27.25<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(10.45,13.37,8.91,9.15,11.94,8.92,11.6,8.78), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")

#Convert blank runs used for background respiration
convert.rMR('empty_pre_9.19.2021_raw.txt', 
            'Converted Empty_Pre 09.19.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_pre_9.21.2021_raw.txt',
            'Converted Empty_Post 09.21.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_pre_9.23.2021_raw.txt',
            'Converted Empty_Post 09.23.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('empty_pre_9.28.2021_raw.txt',
            'Converted Empty_Pre 09.28.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('empty_post_pre_9.30.2021_raw.txt',
            'Converted Empty_Post 09.30.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('empty_post_10.2.2021_raw.txt',
            'Converted Empty_Post 10.02.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('empty_pre_9.24.2021_raw.txt',
            'Converted Empty_Pre 09.24.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")



#Import Background Respiration Files
LF27.pre<-import.test('Converted Empty_Pre 09.19.2021.txt', 
                     info.data=info_27.20, 
                     logger= "AutoResp",
                     n.chamber=8,
                     plot.temperature= F,
                     plot.oxygen = F)
LF27.post<-import.test('Converted Empty_Post 09.21.2021.txt',
                      info.data=info_27.21, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
BC27.pre<-import.test('Converted Empty_Post 09.21.2021.txt',
                      info.data=info_27.21, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
BC27.post<-import.test('Converted Empty_Post 09.23.2021.txt',
              info.data=info_27.23, 
              logger="AutoResp",
              n.chamber=8,
              plot.temperature = F,
              plot.oxygen = F)

DB27.pre<-import.test('Converted Empty_Pre 09.28.2021.txt',
                      info.data=info_27.29, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
DB27.post<-import.test('Converted Empty_Post 09.30.2021.txt',
                       info.data=info_27.30, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)

LFDB31.pre<-import.test('Converted Empty_Post 09.30.2021.txt',
                      info.data=info_31.01, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
LFDB31.post<-import.test('Converted Empty_Post 10.02.2021.txt',
                       info.data=info_31.02, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)

DB227.pre<-import.test('Converted Empty_Pre 09.24.2021.txt',
                      info.data=info_27.29, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
DB227.post<-import.test('Converted Empty_Pre 09.28.2021.txt',
                       info.data=info_27.30, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)


#######################################################################
#Standard Metabolic Rate
#######################################################################
#Convert %Air Saturation to mg/L for SMR measurements and Empty chamber runs
#First Run of Accute SMR Upper Pennisula
convert.rMR('SMR_27_LF_9.20.2021_raw.txt',
            'Converted SMR 27 09.20.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_27_BC_9.22.2021_raw.txt',
            'Converted SMR 27 09.22.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('SMR_27_DB_9.29.2021_raw.txt',
            'Converted SMR 27 09.29.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('SMR_31_LF_DB_10.1.2021_raw.txt',
            'Converted SMR 31 10.01.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('SMR_27_DB_9.25.2021_raw.txt',
            'Converted SMR 27 09.25.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")


#Import Raw SMR Data
LF.SMR.27.raw<-import.meas('Converted SMR 27 09.20.2021.txt', 
                        info.data=info_27.20, 
                        logger="AutoResp", 
                        n.chamber=8,
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)
BC.SMR.27.raw<-import.meas('Converted SMR 27 09.22.2021.txt', 
                          info.data=info_27.22, 
                          logger="AutoResp", 
                          n.chamber=8,
                          date.format="MDY",
                          plot.temperature = F,
                          plot.oxygen = F)
DB.SMR.27.raw<-import.meas('Converted SMR 27 09.29.2021.txt',
                            info.data=info_27.29,
                            logger="AutoResp",
                            n.chamber=8,date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)

SMR.31.raw<-import.meas('Converted SMR 31 10.01.2021.txt',
                           info.data=info_31.01,
                           logger="AutoResp",
                           n.chamber=8,date.format="MDY",
                           plot.temperature = F,
                           plot.oxygen = F)

DB2.SMR.27.raw<-import.meas('Converted SMR 27 09.25.2021.txt',
                           info.data=info_27.25,
                           logger="AutoResp",
                           n.chamber=8,date.format="MDY",
                           plot.temperature = F,
                           plot.oxygen = F)


#Clean Raw SMR data with background respiration measurments
LF.SMR.27.clean<-correct.meas(info.data=info_27.20,
                               pre.data=LF27.pre,
                               post.data=LF27.post,
                           meas.data=LF.SMR.27.raw,
                           method="exponential")
BC.SMR.27.clean<-correct.meas(info.data=info_27.22,
                               pre.data=BC27.pre,
                               post.data=BC27.post, 
                             meas.data=BC.SMR.27.raw,
                             method="exponential")
DB.SMR.27.clean<-correct.meas(info.data=info_27.29,
                               pre.data=DB27.pre,
                               post.data=DB27.post,
                           meas.data=DB.SMR.27.raw,
                           method="exponential")

SMR.31.clean<-correct.meas(info.data=info_31.01,
                               pre.data=LFDB31.pre,
                               post.data=LFDB31.post,
                               meas.data=SMR.31.raw,
                               method="exponential")


DB2.SMR.27.clean<-correct.meas(info.data=info_27.25,
                              pre.data=DB227.pre,
                              post.data=DB227.post,
                              meas.data=DB2.SMR.27.raw,
                              method="exponential")
#Extract slopes for Standard Metabolic Rate
LF.SMR.27.slopes<-extract.slope(LF.SMR.27.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)

BC.SMR.27.slopes<-extract.slope(BC.SMR.27.clean,
                               method="calcSMR.quant",
                               r2=0.95,
                               p=0.2)
DB.SMR.27.slopes<-extract.slope(DB.SMR.27.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)

SMR.31.slopes<-extract.slope(SMR.31.clean,
                                method="calcSMR.quant",
                                r2=0.95,
                                p=0.2)

DB2.SMR.27.slopes<-extract.slope(DB2.SMR.27.clean,
                                method="calcSMR.quant",
                                r2=0.95,
                                p=0.2)


##Extract slopes for Maximum Metabolic Rate during SMR Measurements
LF.SMR.MMR.27.slopes<-extract.slope(LF.SMR.27.clean,
                                 method="max",
                                 n.slope=1)
BC.SMR.MMR.27.slopes<-extract.slope(BC.SMR.27.clean,
                                 method="max",
                                 n.slope=1)
DB.SMR.MMR.27.slopes<-extract.slope(DB.SMR.27.clean,
                                 method="max",
                                 n.slope=1)
SMR.MMR.31.slopes<-extract.slope(SMR.31.clean,
                                    method="max",
                                    n.slope=1)
DB2.SMR.MMR.27.slopes<-extract.slope(DB2.SMR.27.clean,
                                     method="max",
                                     n.slope=1)




##Calculate Standard Metabolic Rates
LF.SMR.27<-calculate.MR(LF.SMR.27.slopes,density=1000)
BC.SMR.27<-calculate.MR(BC.SMR.27.slopes,density=1000)
DB.SMR.27<-calculate.MR(DB.SMR.27.slopes,density=1000)
SMR.31<-calculate.MR(SMR.31.slopes,density=1000)

DB2.SMR.27<-calculate.MR(DB2.SMR.27.slopes,density=1000)

##Calculate Maximum Metabolic Rates 
LF.SMR.MMR.27<-calculate.MR(LF.SMR.MMR.27.slopes,density=1000)
BC.SMR.MMR.27<-calculate.MR(BC.SMR.MMR.27.slopes,density=1000)
DB.SMR.MMR.27<-calculate.MR(DB.SMR.MMR.27.slopes,density=1000)
SMR.MMR.31<-calculate.MR(SMR.MMR.31.slopes,density=1000)

DB2.SMR.MMR.27<-calculate.MR(DB2.SMR.MMR.27.slopes,density=1000)


##Combine Data into one Dataframe
SMR.27.Acute<-rbind.data.frame(LF.SMR.27,BC.SMR.27,DB.SMR.27)
SMR.MMR.27.Acute<-rbind.data.frame(LF.SMR.MMR.27,BC.SMR.MMR.27,DB.SMR.MMR.27)

DB2<-cbind(DB2.SMR.27,Trial,Temperature,Event)
DB2.SMR.MMR<-cbind(DB2.SMR.MMR.27,Trial,Temperature,Event)

write.csv(DB2, "DB2 27 Acute SMR.csv")
write.csv(DB2.SMR.MMR, "DB2 27 Acute MMR.SMR.csv")

##Add additional Columns
Temperature<-27
Trial<-"SMR"
Event<-1

SMR.27.Acute<-cbind(SMR.27.Acute,Trial,Temperature,Event)
SMR.MMR.27.Acute<-cbind(SMR.MMR.27.Acute,Trial,Temperature,Event)
LF.27.SMR.MMR.Acute<-cbind(SMR.MMR.27.Acute,Trial,Temperature,Event)
LF.27.SMR.Acute<-cbind(LF.SMR.27,Trial,Temperature,Event)
Temperature<-31
SMR.MMR.31.Acute<-cbind(SMR.MMR.31,Trial,Temperature,Event)
SMR.31.Acute<-cbind(SMR.31,Trial,Temperature, Event)
##############33
##Save out CSVs

write.csv(SMR.27.Acute, "All 27 Acute SMR.csv")
write.csv(SMR.MMR.27.Acute, "All 27 Acute SMR.MMR.csv")
write.csv(LF.27.SMR.Acute,"LF 27 Acute SMR.csv")
write.csv(SMR.MMR.31.Acute, "31 Acute SMR-MMR.csv")
write.csv(SMR.31.Acute,"SMR 31 Acute.csv")
write.csv(LF.27.SMR.MMR.Acute,"LF 27 Acute SMR-MMR.csv")


#######################################################################
#Maximum Metabolic Rate
########################################################################

#Convert %Air Saturation to mg/L for MMR measurements and Empty chamber runs
#First Run of Acute MMR Upper Pennisula
convert.rMR('MMR_27_LF_9.21.2021_raw.txt',
            'Converted MMR 27 09.21.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_27_BC_9.23.2021_raw.txt',
            'Converted MMR 27 09.23.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_27_DB_9.30.2021_raw.txt',
            'Converted MMR 27 09.30.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('MMR_31_LF_DB_10.2.2021_raw.txt',
            'Converted MMR 31 10.02.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw MMR Data
LF.MMR.27.raw<-import.meas('Converted MMR 27 09.21.2021.txt', 
                            info.data=info_27.21, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
BC.MMR.27.raw<-import.meas('Converted MMR 27 09.23.2021.txt', 
                            info.data=info_27.23, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
DB.MMR.27.raw<-import.meas('Converted MMR 27 09.30.2021.txt', 
                            info.data=info_27.30, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
MMR.31.raw<-import.meas('Converted MMR 31 10.02.2021.txt', 
                            info.data=info_27.02, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)


#Clean Raw MMR data with background respiration measurments
LF.MMR.27.clean<-correct.meas(info.data=info_27.21,
                               pre.data=LF27.pre,
                               post.data=LF27.post,
                               meas.data=LF.MMR.27.raw,
                               method="exponential")
                                
BC.MMR.27.clean<-correct.meas(info.data=info_27.23,
                               pre.data=BC27.pre,
                               post.data=BC27.post, 
                               meas.data=BC.MMR.27.raw,
                               method="exponential")
DB.MMR.27.clean<-correct.meas(info.data=info_27.30,
                               pre.data=DB27.pre,
                               post.data=DB27.post,
                               meas.data=DB.MMR.27.raw,
                               method="exponential")
MMR.31.clean<-correct.meas(info.data=info_31.02,
                              pre.data=LFDB31.pre,
                              post.data=LFDB31.post,
                              meas.data=MMR.31.raw,
                              method="exponential")

#Extract slopes
LF.MMR.27.slopes<-extract.slope(LF.MMR.27.clean,
                                 method="max",
                                 n.slope=1)
BC.MMR.27.slopes<-extract.slope(BC.MMR.27.clean,
                                 method="max",
                                 n.slope=1)
DB.MMR.27.slopes<-extract.slope(DB.MMR.27.clean,
                                 method="max",
                                 n.slope=1)
MMR.31.slopes<-extract.slope(MMR.31.clean,
                                method="max",
                                n.slope=1)

##Calculate Metabolic Rates
LF.MMR.27<-calculate.MR(LF.MMR.27.slopes,density=1000)
BC.MMR.27<-calculate.MR(BC.MMR.27.slopes,density=1000)
DB.MMR.27<-calculate.MR(DB.MMR.27.slopes,density=1000)
MMR.31<-calculate.MR(MMR.31.slopes,density=1000)


##Combine data into one Dataframe
LF.27.MMR.Acute<-rbind.data.frame(LF.MMR.27)
MMR.27.Acute<-rbind.data.frame(LF.MMR.27,BC.MMR.27,DB.MMR.27)
##Add additional columns
Trial<-"MMR"
LF.27.MMR.Acute<-cbind(LF.27.MMR.Acute,Trial,Temperature,Event)
MMR.27.Acute$Trial<-"MMR"
MMR.27.Acute$Temperature<-27
MMR.27.Acute$Event<-1
MMR.31$Trial<-"MMR"
MMR.31$Temperature<-31
MMR.31$Event<-1
##Save data our as CSV
write.csv(MMR.27.Acute,"All 27 Acute MMR.csv")
write.csv(MMR.31,"31 Acute MMR.csv")
write.csv(LF.27.MMR.Acute,"LF 27 Acute MMR.csv")

#######################################################################
#Data Analysis and grapahs
#######################################################################


All.27.MMR<-read.csv("All 27 Acute MMR.csv")
All.27.SMR<-read.csv("All 27 Acute SMR.csv")
Best<-read.csv("27 Acute MMR Best.csv")
SMR.31.Acute

s<-ggplot(DB2.SMR,aes(Trial, MR.mass/60, group=Ind,color=Ind)) 
s<-s+geom_line() + geom_point() + ggtitle("DB2 SMR, Comparison")+theme_classic()
s<-s+ylab("MO2 (mgO2/(kg min))")
s
DB2.SMR.Acute$Temperature<-factor(Temperature)

SMR<-ggplot(DB2.SMR.Acute,aes(x=Temperature,y=MR.mass/60))
SMR<-SMR+geom_boxplot() + ggtitle("DB4 SMR, SMR")+theme_classic()
SMR<-SMR+geom_boxplot(data=DB2.MMR, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ylab("MO2 (mgO2/(kg min))")
SMR


grid.arrange(
ggplot(All.27.SMR,aes(as.factor(Temperature),y=MR.mass/60)) + 
                geom_boxplot() +geom_boxplot(data=All.27.MMR, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
                ylab("MO2 (mgO2/(kg min))") + ggtitle("DB SMR,Best MMR ") + theme(legend.position="none"),
ggplot(All.27.SMR,aes(x=as.factor(Temperature),y=MR.mass/60)) + 
        geom_boxplot() +geom_boxplot(data=Best, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
        ylab("MO2 (mgO2/(kg min))") + ggtitle("DB SMR, MMR ") + theme(legend.position="none"),
ggplot(SMR.MMR.27.Acute,aes(x=as.factor(Temperature),y=MR.mass/60)) + 
        geom_boxplot() +geom_boxplot(data=SMR.MMR.27.Acute, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
        ylab("MO2 (mgO2/(kg min))") + ggtitle("DB SMR, SMR-MMR ") + theme(legend.position="none"),
        
        ncol=3
)

grid.arrange(
        ggplot(subset(All.27.MMR,Event%in%c(1)), 
               aes(x = Trial, y = MR.mass, group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, DB 25") + theme(legend.position="none")+ylim(200,1500),
        ggplot(subset(DB2.MMR,Event%in%c(4)), 
               aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, DB 25.2") + theme(legend.position="none")+ylim(200,1500),
        ggplot(subset(DB2.MMR,Event%in%c(2)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, DB 23") + theme(legend.position="none")+ylim(200,1500),
        ggplot(subset(DB2.MMR,Event%in%c(3)), 
               aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, DB 21") + theme(legend.position="none")+ylim(200,1500),
        
        ncol=2
)
aMMR<-aov(MMRC$MR.mass~MMRC$Temperature*MMRC$Trial,data=MMRC)
Anova(aMMR)
TukeyHSD((aMMR))

as.factor(DB2.MMR$Temperature)


