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
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_9_6_2021")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used

info_21.06<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(4.37,4.91,5.24,4.86,4.92,4.62,4.78,3.76), 
                       Volume = c(78,78,78,78,78,78,78,78), 
                       DO.unit = "mg/L")
info_21.07<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(4.17,4.93,4.92,4.68,4.64,4.48,4.57,3.59), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_23.08<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(4.22,4.64,4.84,4.55,4.62,4.42,4.52,3.55), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_25.09<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                                Mass=c(4.02,4.56,4.75,4.51,4.36,4.38,4.41,3.48), 
                                Volume = c(78,78,78,78,78,78,78,78), 
                                DO.unit = "mg/L")
info_21.2.10<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                        Mass=c(4.09,4.54,4.7,4.47,4.36,4.38,4.45,3.41), 
                        Volume = c(78,78,78,78,78,78,78,78), 
                        DO.unit = "mg/L")

#Convert blank runs used for background respiration
convert.rMR('empty_pre_9.6.2021_raw.txt', 
            'Converted Empty_Pre 09.06.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_9.10.2021_raw.txt',
            'Converted Empty_Post 09.10.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")


#Import Background Respiration Files
DB1.pre<-import.test('Converted Empty_Pre 09.06.2021.txt', 
                     info.data=info_21.04, 
                     logger= "AutoResp",
                     n.chamber=8,
                     plot.temperature= F,
                     plot.oxygen = F)
DB1.post<-import.test('Converted Empty_Post 09.10.2021.txt',
                      info.data=info_21.04, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)

#######################################################################
#Standard Metabolic Rate
#######################################################################
#Convert %Air Saturation to mg/L for SMR measurements and Empty chamber runs
#First Run of Accute SMR Upper Pennisula
convert.rMR('SMR_21_DB_9.6.2021_raw.txt',
            'Converted SMR 21 09.06.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_21.2_DB_9.9.2021_raw.txt',
            'Converted SMR 21.2 09.09.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_23_DB_9.7.2021_raw.txt',
            'Converted SMR 23 09.07.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_25_DB_9.8.2021_raw.txt',
            'Converted SMR 25 09.08.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw SMR Data
DB1.SMR.21.raw<-import.meas('Converted SMR 21 09.06.2021.txt', 
                        info.data=info_21.06, 
                        logger="AutoResp", 
                        n.chamber=8,
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)
DB1.SMR.21.2.raw<-import.meas('Converted SMR 21.2 09.09.2021.txt', 
                        info.data=info_25.08, 
                        n.chamber=8,
                        logger="AutoResp", 
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)
DB1.SMR.23.raw<-import.meas('Converted SMR 23 09.07.2021.txt', 
                          info.data=info_21.07, 
                          logger="AutoResp", 
                          n.chamber=8,
                          date.format="MDY",
                          plot.temperature = F,
                          plot.oxygen = F)
DB1.SMR.25.raw<-import.meas('Converted SMR 25 09.08.2021.txt', 
                        info.data=info_23.08, 
                        logger="AutoResp", 
                        n.chamber=8,
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)

#Clean Raw SMR data with background respiration measurments
DB1.SMR.21.clean<-correct.meas(info.data=info_21.06,
                               pre.data=DB1.pre,
                               post.data=DB1.post,
                           meas.data=DB1.SMR.21.raw,
                           method="exponential")
DB1.SMR.21.2.clean<-correct.meas(info.data=info_25.09,
                           pre.data=DB1.pre,
                           post.data=DB1.post, 
                           meas.data=DB1.SMR.21.2.raw,
                           method="exponential")
DB1.SMR.23.clean<-correct.meas(info.data=info_21.07,
                               pre.data=DB1.pre,
                               post.data=DB1.post, 
                             meas.data=DB1.SMR.23.raw,
                             method="exponential")
DB1.SMR.25.clean<-correct.meas(info.data=info_23.08,
                               pre.data=DB1.pre,
                               post.data=DB1.post,
                           meas.data=DB1.SMR.25.raw,
                           method="exponential")

##Extract slopes for Standard Metabolic Rate
DB1.SMR.21.slopes<-extract.slope(DB1.SMR.21.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)
DB1.SMR.21.2.slopes<-extract.slope(DB1.SMR.21.2.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)
DB1.SMR.23.slopes<-extract.slope(DB1.SMR.23.clean,
                               method="calcSMR.quant",
                               r2=0.95,
                               p=0.2)
DB1.SMR.25.slopes<-extract.slope(DB1.SMR.25.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)
##Extract Slopes for Maximum Metabolic Rates from SMR Data
DB1.SMR.MMR.21.slopes<-extract.slope(DB1.SMR.21.clean,
                                     method="max",
                                     n.slope=1)

DB1.SMR.MMR.21.2.slopes<-extract.slope(DB1.SMR.21.2.clean,
                                       method="max",
                                       n.slope=1)

DB1.SMR.MMR.23.slopes<-extract.slope(DB1.SMR.23.clean,
                                     method="max",
                                     n.slope=1)

DB1.SMR.MMR.25.slopes<-extract.slope(DB1.SMR.25.clean,
                                     method="max",
                                     n.slope=1)

##Calculate Standard Metabolic Rates
DB1.SMR.21<-calculate.MR(DB1.SMR.21.slopes,density=1000)
DB1.SMR.21.2<-calculate.MR(DB1.SMR.21.2.slopes,density=1000)
DB1.SMR.23<-calculate.MR(DB1.SMR.23.slopes,density=1000)
DB1.SMR.25<-calculate.MR(DB1.SMR.25.slopes,density=1000)


##Calculate Maximum Metabolic Rates
DB1.SMR.MMR.21<-calculate.MR(DB1.SMR.MMR.21.slopes,density=1000)
DB1.SMR.MMR.21.2<-calculate.MR(DB1.SMR.MMR.21.2.slopes,density=1000)
DB1.SMR.MMR.23<-calculate.MR(DB1.SMR.MMR.23.slopes,density=1000)
DB1.SMR.MMR.25<-calculate.MR(DB1.SMR.MMR.25.slopes,density=1000)


##Combine Data into one dataframe
DB1.SMR.Acute<-rbind.data.frame(DB1.SMR.21,DB1.SMR.23,DB1.SMR.25,DB1.SMR.21.2)
DB1.SMR.MMR.Acute<-rbind.data.frame(DB1.SMR.MMR.21,DB1.SMR.MMR.23,DB1.SMR.MMR.25,DB1.SMR.MMR.21.2)

##Add additional columns
#Temperature as a Factor
Temperature<-c(21,21,21,21,21,21,21,21,23,23,23,23,23,23,23,23,25,25,25,25,25,25,25,25,21,21,21,21,21,21,21,21)
#Trial measurement for SMR.MMR
Trial<-c("SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR","SMR","SMR")
#Experiment Order
Event<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4)

DB1.SMR.MMR.Acute<-cbind(DB1.SMR.MMR.Acute,Trial,Temperature,Event)
DB1.SMR.Acute<-cbind(DB1.SMR.Acute,Temperature, Event)

########################################################################
##Save out CSVs
setwd("~/Desktop/Dissertation/Chapter 2")

write.csv(DB1.SMR.Acute,"LP DB Acute SMR EXP 1.csv")
write.csv(DB1.SMR.MMR.Acute,"LP DB Acute SMR.MMR EXP 1.csv")

#######################################################################
#Maximum Metabolic Rate
########################################################################

#Convert %Air Saturation to mg/L for MMR measurements and Empty chamber runs
#First Run of Acute MMR Upper Pennisula
convert.rMR('MMR_21_DB_9.7.2021_raw.txt',
            'Converted MMR 21 09.07.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_21.2_DB_9.10.2021_raw.txt',
            'Converted MMR 21.2 09.10.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_23_DB_9.8.2021_raw.txt',
            'Converted MMR 23 09.08.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_25_DB_9.9.2021_raw.txt',
            'Converted MMR 25 09.09.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw MMR Data
DB1.MMR.21.raw<-import.meas('Converted MMR 21 09.07.2021.txt', 
                            info.data=info_21.07, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
DB1.MMR.21.2.raw<-import.meas('Converted MMR 21.2 09.10.2021.txt', 
                              info.data=info_21.2.10, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
DB1.MMR.23.raw<-import.meas('Converted MMR 23 09.08.2021.txt', 
                            info.data=info_23.08, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
DB1.MMR.25.raw<-import.meas('Converted MMR 25 09.09.2021.txt', 
                            info.data=info_25.09, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)

#Clean Raw MMR data with background respiration measurments
DB1.MMR.21.clean<-correct.meas(info.data=info_21.07,
                               pre.data=DB1.pre,
                               post.data=DB1.post,
                               meas.data=DB1.MMR.21.raw,
                               method="exponential")
DB1.MMR.21.2.clean<-correct.meas(info.data=info_21.2.10,
                                 pre.data=DB1.pre,
                                 post.data=DB1.post, 
                                 meas.data=DB1.MMR.21.2.raw,
                                 method="exponential")
DB1.MMR.23.clean<-correct.meas(info.data=info_23.08,
                               pre.data=DB1.pre,
                               post.data=DB1.post, 
                               meas.data=DB1.MMR.23.raw,
                               method="exponential")
DB1.MMR.25.clean<-correct.meas(info.data=info_25.09,
                               pre.data=DB1.pre,
                               post.data=DB1.post,
                               meas.data=DB1.MMR.25.raw,
                               method="exponential")

#Extract slopes
DB1.MMR.21.slopes<-extract.slope(DB1.MMR.21.clean,
                                 method="max",
                                 n.slope=1)
DB1.MMR.21.2.slopes<-extract.slope(DB1.MMR.21.2.clean,
                                   method="max",
                                   n.slope=1)
DB1.MMR.23.slopes<-extract.slope(DB1.MMR.23.clean,
                                 method="max",
                                 n.slope=1)
DB1.MMR.25.slopes<-extract.slope(DB1.MMR.25.clean,
                                 method="max",
                                 n.slope=1)

##Calculate Metabolic Rates
DB1.MMR.21<-calculate.MR(DB1.MMR.21.slopes,density=1000)
DB1.MMR.21.2<-calculate.MR(DB1.MMR.21.2.slopes,density=1000)
DB1.MMR.23<-calculate.MR(DB1.MMR.23.slopes,density=1000)
DB1.MMR.25<-calculate.MR(DB1.MMR.25.slopes,density=1000)


##Combine into one Datafram
DB1.MMR.Acute<-rbind.data.frame(DB1.MMR.21,DB1.MMR.23,DB1.MMR.25,DB1.MMR.21.2)

##Add columns
Trial<-c("MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR","MMR","MMR")
DB1.MMR.Acute<-cbind(DB1.MMR.Acute,Trial,Temperature,Event)

#################3
##Save out as CSV

setwd("~/Desktop/Dissertation/Chapter 2")

write.csv(DB1.MMR.Acute,"LP DB Acute MMR EXP 1.csv")

#######################################################################
#Data Analysis and grapahs
#######################################################################

DB1.MMR<-rbind(DB1.SMR.MMR.Acute,DB1.MMR.Acute)
DB1.MMR<-read.csv("LP DB MMR BEST EXP 1.csv")
DB1.MMR$Temperature<-factor(Temperature)
DB1.SMR.MMR.Acute$Temperature<-factor(Temperature)
DB1.MMR.Acute$Temperature<-factor(Temperature)
DB1.SMR.Acute$Temperature<-factor(Temperature)
MMRC$Temperature<-factor(Temperature)
MMRC<-as.data.frame(MMRC)
m<-ggplot() 
m<-m+geom_point(data=DB1.SMR.MMR.Acute,aes(Temperature, MR.mass/60, group=Ind,color=Ind)) + ggtitle("DB1 MMR, Comparison")+theme_classic()
m<-m+geom_point(data=DB1.MMR.Acute,aes(Temperature, MR.mass/60, group=Ind,color=Ind))
m<-m+ylab("MO2 (mgO2/(kg min))")
m

s<-ggplot(DB2.SMR,aes(Trial, MR.mass/60, group=Ind,color=Ind)) 
s<-s+geom_line() + geom_point() + ggtitle("DB2 SMR, Comparison")+theme_classic()
s<-s+ylab("MO2 (mgO2/(kg min))")




DB1.SMR.Acute$Stock <-  "LP"
LF3.SMR.Acute$Stock<-"UP"

SMR.Acute<-rbind(DB1.SMR.Acute,LF3.SMR.Acute)

DB1.MMR$Stock <- "LP"
LF3.MMR$Stock<-"UP"
DB1.MMR.Acute$Stock<- "LP"
LF3.MMR.Acute$Stock<-"UP"
MMR.Acute<-rbind(DB1.MMR,LF3.MMR)
DB1.SMR.MMR.Acute$Stock<-"LP"
MMR<-rbind(DB1.MMR.Acute,DB1.SMR.MMR.Acute)


grid.arrange(
        ggplot(DB1.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
                geom_boxplot() +geom_boxplot(data=DB1.MMR, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
                ylab("MO2 (mgO2/(kg min))") + ggtitle("DB SMR,Best MMR ") + theme(legend.position="none"),
        ggplot(DB1.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
                geom_boxplot() +geom_boxplot(data=DB1.MMR.Acute, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
                ylab("MO2 (mgO2/(kg min))") + ggtitle("DB SMR, MMR ") + theme(legend.position="none"),
        ggplot(DB1.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
                geom_boxplot() +geom_boxplot(data=DB1.SMR.MMR.Acute, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
                ylab("MO2 (mgO2/(kg min))") + ggtitle("DB SMR, SMR-MMR ") + theme(legend.position="none"),
        
        ncol=3
)





grid.arrange(
        ggplot(subset(DB1.MMR,Event%in%c(1)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, DB 21") + theme(legend.position="none"),
        ggplot(subset(DB1.MMR,Event%in%c(4)), 
               aes(x = Trial, y = MR.mass,  group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, DB 21.2") + theme(legend.position="none"),
        ggplot(subset(DB1.MMR,Event%in%c(2)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, DB 23") + theme(legend.position="none"),
        ggplot(subset(DB1.MMR,Event%in%c(3)), 
               aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, DB 25") + theme(legend.position="none"),
        ncol=4
        
        
)





