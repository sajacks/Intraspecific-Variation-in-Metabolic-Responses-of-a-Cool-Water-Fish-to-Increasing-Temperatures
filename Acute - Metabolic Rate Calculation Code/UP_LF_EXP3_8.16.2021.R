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
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_8.16.2021")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used
info_21.16<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(2.59,4.42,4.63,4.39,4.28,4.88,4.24,4.97), 
                       Volume = c(78,78,78,78,78,78,78,78), 
                       DO.unit = "mg/L")
info_21.17<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(2.44,4.34,4.38,4.1,4.25,4.64,4.15,4.81), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_23.18<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(2.42, 4.21,4.41,4.01,4.05,4.53,4.05,4.73), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_25.19<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                                Mass=c(2.35,4.18,4.28,3.98,4.12,4.53,4.04,4.63), 
                                Volume = c(78,78,78,78,78,78,78,78), 
                                DO.unit = "mg/L")
info_21.2.20<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                        Mass=c(2.39,4.07,4.26,3.95,3.99,4.55,4.08,4.67), 
                        Volume = c(78,78,78,78,78,78,78,78), 
                        DO.unit = "mg/L")

#Convert blank runs used for background respiration
convert.rMR('empty_pre_8.16.2021_raw.txt', 
            'Converted Empty_Pre 08.16.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_8.20.2021_raw.txt',
            'Converted Empty_Post 08.20.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")


#Import Background Respiration Files
LF3.pre<-import.test('Converted Empty_Pre 08.16.2021.txt', 
                     info.data=info_21.16, 
                     logger= "AutoResp",
                     n.chamber=8,
                     plot.temperature= F,
                     plot.oxygen = F)
LF3.post<-import.test('Converted Empty_Post 08.20.2021.txt',
                      info.data=info_21.16, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)

#######################################################################
#Standard Metabolic Rate
#######################################################################
#Convert %Air Saturation to mg/L for SMR measurements and Empty chamber runs
#First Run of Accute SMR Upper Pennisula
convert.rMR('SMR_21_LF_8.16.2021_raw.txt',
            'Converted SMR 21 08.16.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_21.2_LF_8.19.2021_raw.txt',
            'Converted SMR 21.2 08.19.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_23_LF_8.17.2021_raw.txt',
            'Converted SMR 23 08.17.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_25_LF_8.18.2021_raw.txt',
            'Converted SMR 25 08.18.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw SMR Data
LF3.SMR.21.raw<-import.meas('Converted SMR 21 08.16.2021.txt', 
                        info.data=info_21.16, 
                        logger="AutoResp", 
                        n.chamber=8,
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)
LF3.SMR.21.2.raw<-import.meas('Converted SMR 21.2 08.19.2021.txt', 
                        info.data=info_25.19, 
                        n.chamber=8,
                        logger="AutoResp", 
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)
LF3.SMR.23.raw<-import.meas('Converted SMR 23 08.17.2021.txt', 
                          info.data=info_21.17, 
                          logger="AutoResp", 
                          n.chamber=8,
                          date.format="MDY",
                          plot.temperature = F,
                          plot.oxygen = F)
LF3.SMR.25.raw<-import.meas('Converted SMR 25 08.18.2021.txt',
                            info.data=info_23.18,
                            logger="AutoResp",
                            n.chamber=8,date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)

#Clean Raw SMR data with background respiration measurments
LF3.SMR.21.clean<-correct.meas(info.data=info_21.16,
                               pre.data=LF3.pre,
                               post.data=LF3.post,
                           meas.data=LF3.SMR.21.raw,
                           method="exponential")
LF3.SMR.21.2.clean<-correct.meas(info.data=info_21.2.19,
                           pre.data=LF3.pre,
                           post.data=LF3.post, 
                           meas.data=LF3.SMR.21.2.raw,
                           method="exponential")
LF3.SMR.23.clean<-correct.meas(info.data=info_23.17,
                               pre.data=LF3.pre,
                               post.data=LF3.post, 
                             meas.data=LF3.SMR.23.raw,
                             method="exponential")
LF3.SMR.25.clean<-correct.meas(info.data=info_25.18,
                               pre.data=LF3.pre,
                               post.data=LF3.post,
                           meas.data=LF3.SMR.25.raw,
                           method="exponential")

#Extract slopes for Standard Metabolic Rate
LF3.SMR.21.slopes<-extract.slope(LF3.SMR.21.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)
LF3.SMR.21.2.slopes<-extract.slope(LF3.SMR.21.2.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)
LF3.SMR.23.slopes<-extract.slope(LF3.SMR.23.clean,
                               method="calcSMR.quant",
                               r2=0.95,
                               p=0.2)
LF3.SMR.25.slopes<-extract.slope(LF3.SMR.25.clean,
                             method="calcSMR.quant",
                             r2=0.9,
                             p=0.2)

##Extract slopes for Maximum Metabolic Rate during SMR Measurements
LF3.SMR.MMR.21.slopes<-extract.slope(LF3.SMR.21.clean,
                                 method="max",
                                 n.slope=1)
LF3.SMR.MMR.21.2.slopes<-extract.slope(LF3.SMR.21.2.clean,
                                   method="max",
                                   n.slope=1)
LF3.SMR.MMR.23.slopes<-extract.slope(LF3.SMR.23.clean,
                                 method="max",
                                 n.slope=1)
LF3.SMR.MMR.25.slopes<-extract.slope(LF3.SMR.25.clean,
                                 method="max",
                                 n.slope=1)



##Calculate Standard Metabolic Rates
LF3.SMR.21<-calculate.MR(LF3.SMR.21.slopes,density=1000)
LF3.SMR.21.2<-calculate.MR(LF3.SMR.21.2.slopes,density=1000)
LF3.SMR.23<-calculate.MR(LF3.SMR.23.slopes,density=1000)
LF3.SMR.25<-calculate.MR(LF3.SMR.25.slopes,density=1000)

##Calculate Maximum Metabolic Rates 
LF3.SMR.MMR.21<-calculate.MR(LF3.SMR.MMR.21.slopes,density=1000)
LF3.SMR.MMR.21.2<-calculate.MR(LF3.SMR.MMR.21.2.slopes,density=1000)
LF3.SMR.MMR.23<-calculate.MR(LF3.SMR.MMR.23.slopes,density=1000)
LF3.SMR.MMR.25<-calculate.MR(LF3.SMR.MMR.25.slopes,density=1000)

##Combine Data into one Dataframe
LF3.SMR.Acute<-rbind.data.frame(LF3.SMR.21,LF3.SMR.21.2,LF3.SMR.23,LF3.SMR.25)
LF3.SMR.MMR.Acute<-rbind.data.frame(LF3.SMR.MMR.21,LF3.SMR.MMR.21.2,LF3.SMR.MMR.23,LF3.SMR.MMR.25)

##Add additional Columns
Temperature<-c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,23,23,23,23,23,23,23,23,25,25,25,25,25,25,25,25)
Trial<-c("SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR","SMR","SMR")
Event<-c(1,1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3)

LF3.SMR.Acute<-cbind(LF3.SMR.Acute,Trial,Temperature,Event)
LF3.SMR.MMR.Acute<-cbind(LF3.SMR.MMR.Acute,Trial,Temperature,Event)

##############33
##Save out CSVs

setwd("~/Desktop/Dissertation/Chapter 2")
write.csv(LF3.SMR.Acute,"UP LF Acute SMR EXP 3.csv")
write.csv(LF3.SMR.MMR.Acute,"UP LF Acute SMR-MMR EXP 3.csv")


#######################################################################
#Maximum Metabolic Rate
########################################################################

#Convert %Air Saturation to mg/L for MMR measurements and Empty chamber runs
#First Run of Acute MMR Upper Pennisula
convert.rMR('MMR_21_LF_8.17.2021_raw.txt',
            'Converted MMR 21 08.17.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_21.2_LF_8.20.2021_raw.txt',
            'Converted MMR 21.2 08.20.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_23_LF_8.18.2021_raw.txt',
            'Converted MMR 23 08.18.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_25_LF_8.19.2021_raw.txt',
            'Converted MMR 25 08.19.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw MMR Data
LF3.MMR.21.raw<-import.meas('Converted MMR 21 08.17.2021.txt', 
                            info.data=info_21.16, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
LF3.MMR.21.2.raw<-import.meas('Converted MMR 21.2 08.20.2021.txt', 
                              info.data=info_21.2.19, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
LF3.MMR.23.raw<-import.meas('Converted MMR 23 08.18.2021.txt', 
                            info.data=info_23.17, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
LF3.MMR.25.raw<-import.meas('Converted MMR 25 08.19.2021.txt', 
                            info.data=info_25.18, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)

#Clean Raw MMR data with background respiration measurments
LF3.MMR.21.clean<-correct.meas(info.data=info_21.16,
                               pre.data=LF3.pre,
                               post.data=LF3.post,
                               meas.data=LF3.MMR.21.raw,
                               method="exponential")
LF3.MMR.21.2.clean<-correct.meas(info.data=info_21.2.19,
                                 pre.data=LF3.pre,
                                 post.data=LF3.post, 
                                 meas.data=LF3.MMR.21.2.raw,
                                 method="exponential")
LF3.MMR.23.clean<-correct.meas(info.data=info_23.17,
                               pre.data=LF3.pre,
                               post.data=LF3.post, 
                               meas.data=LF3.MMR.23.raw,
                               method="exponential")
LF3.MMR.25.clean<-correct.meas(info.data=info_25.18,
                               pre.data=LF3.pre,
                               post.data=LF3.post,
                               meas.data=LF3.MMR.25.raw,
                               method="exponential")

#Extract slopes
LF3.MMR.21.slopes<-extract.slope(LF3.MMR.21.clean,
                                 method="max",
                                 n.slope=1)
LF3.MMR.21.2.slopes<-extract.slope(LF3.MMR.21.2.clean,
                                   method="max",
                                   n.slope=1)
LF3.MMR.23.slopes<-extract.slope(LF3.MMR.23.clean,
                                 method="max",
                                 n.slope=1)
LF3.MMR.25.slopes<-extract.slope(LF3.MMR.25.clean,
                                 method="max",
                                 n.slope=1)

##Calculate Metabolic Rates
LF3.MMR.21<-calculate.MR(LF3.MMR.21.slopes,density=1000)
LF3.MMR.21.2<-calculate.MR(LF3.MMR.21.2.slopes,density=1000)
LF3.MMR.23<-calculate.MR(LF3.MMR.23.slopes,density=1000)
LF3.MMR.25<-calculate.MR(LF3.MMR.25.slopes,density=1000)


##Combine data into one Dataframe
LF3.MMR.Acute<-rbind.data.frame(LF3.MMR.21,LF3.MMR.21.2,LF3.MMR.23,LF3.MMR.25)

##Add additional columns
Trial<-c("MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR","MMR","MMR")
LF3.MMR.Acute<-cbind(LF3.MMR.Acute,Trial,Temperature,Event)

##Save data our as CSV
setwd("~/Desktop/Dissertation/Chapter 2")

write.csv(LF3.MMR.Acute,"UP LF Acute MMR EXP 3.csv")

#######################################################################
#Data Analysis and grapahs
#######################################################################

LF3.MMR<-rbind(LF3.SMR.MMR.Acute,LF3.MMR.Acute)

LF3.MMR<-read.csv("UP LF Acute MMR BEST EXP 3.csv")
LF3.MMR$Temperature<-factor(Temperature)
LF3.SMR.MMR.Acute$Temperature<-factor(Temperature)
LF3.MMR.Acute$Temperature<-factor(Temperature)
MMRC$Temperature<-factor(Temperature)
MMRC<-as.data.frame(MMRC)
m<-ggplot() 
m<-m+geom_point(data=LF3.SMR.MMR.Acute,aes(Temperature, MR.mass/60, group=Ind,color=Ind)) + ggtitle("LF3 MMR, Comparison")+theme_classic()
m<-m+geom_point(data=LF3.MMR.Acute,aes(Temperature, MR.mass/60, group=Ind,color=Ind))
m<-m+ylab("MO2 (mgO2/(kg min))")
m

s<-ggplot(LF2.SMR,aes(Trial, MR.mass/60, group=Ind,color=Ind)) 
s<-s+geom_line() + geom_point() + ggtitle("LF2 SMR, Comparison")+theme_classic()
s<-s+ylab("MO2 (mgO2/(kg min))")
s
LF3.SMR.Acute$Temperature<-factor(Temperature)

SMR<-ggplot(LF3.SMR.Acute,aes(x=Temperature,y=MR.mass/60))
SMR<-SMR+geom_boxplot() + ggtitle("LF3 SMR, SMR")+theme_classic()
SMR<-SMR+geom_boxplot(data=LF3.MMR, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ylab("MO2 (mgO2/(kg min))")
SMR


grid.arrange(
        ggplot(LF3.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
                geom_boxplot() +geom_boxplot(data=LF3.MMR, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
                ylab("MO2 (mgO2/(kg min))") + ggtitle("LF SMR,Best MMR ") + theme(legend.position="none"),
ggplot(LF3.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
        geom_boxplot() +geom_boxplot(data=LF3.MMR.Acute, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
        ylab("MO2 (mgO2/(kg min))") + ggtitle("LF SMR, MMR ") + theme(legend.position="none"),
ggplot(LF3.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
        geom_boxplot() +geom_boxplot(data=LF3.SMR.MMR.Acute, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
        ylab("MO2 (mgO2/(kg min))") + ggtitle("LF SMR, SMR-MMR ") + theme(legend.position="none"),
        
        ncol=3
)

grid.arrange(
        ggplot(subset(MMRC,Day%in%c(1)), 
               aes(x = Trial, y = MR.mass, group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF 21") + theme(legend.position="none"),
        ggplot(subset(MMRC,Day%in%c(4)), 
               aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF 21.2") + theme(legend.position="none"),
        ggplot(subset(MMRC,Day%in%c(2)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF 23") + theme(legend.position="none"),
        ggplot(subset(MMRC,Day%in%c(3)), 
               aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF 25") + theme(legend.position="none"),
        
        ncol=2
)
aMMR<-aov(MMRC$MR.mass~MMRC$Temperature*MMRC$Trial,data=MMRC)
Anova(aMMR)
TukeyHSD((aMMR))

as.factor(LF3.MMR$Temperature)


