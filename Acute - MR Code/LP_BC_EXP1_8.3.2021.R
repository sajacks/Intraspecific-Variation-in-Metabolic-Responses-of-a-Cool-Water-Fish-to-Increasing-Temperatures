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

#Convert blank runs used for background respiration
convert.rMR('empty_pre_8.3.2021_raw.txt', 
            'Converted Empty_Pre 08.03.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_8.7.2021_raw.txt',
            'Converted Empty_Post 08.07.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")


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

#######################################################################
#Standard Metabolic Rate
#######################################################################
#Convert %Air Saturation to mg/L for SMR measurements and Empty chamber runs
#First Run of Accute SMR Upper Pennisula
convert.rMR('SMR_21_BC_8.3.2021_raw.txt',
            'Converted SMR 21 08.03.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR2_21.2_BC_8.6.2021_raw.txt',
            'Converted SMR 21.2 08.06.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_23_BC_8.4.2021_raw.txt',
            'Converted SMR 23 08.04.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMRC_25_BC_8.5.2021_raw.txt',
            'Converted SMR 25 08.05.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

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

##Extract slopes for Standard Metabolic Rate
BC1.SMR.21.slopes<-extract.slope(BC1.SMR.21.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)
BC1.SMR.21.2.slopes<-extract.slope(BC1.SMR.21.2.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)
BC1.SMR.23.slopes<-extract.slope(BC1.SMR.23.clean,
                               method="calcSMR.quant",
                               r2=0.95,
                               p=0.2)
BC1.SMR.25.slopes<-extract.slope(BC1.SMR.25.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)
##Extract Slopes for Maximum Metabolic Rates from SMR Data


##Calculate Standard Metabolic Rates
BC1.SMR.21<-calculate.MR(BC1.SMR.21.slopes,density=1000)
BC1.SMR.21.2<-calculate.MR(BC1.SMR.21.2.slopes,density=1000)
BC1.SMR.23<-calculate.MR(BC1.SMR.23.slopes,density=1000)
BC1.SMR.25<-calculate.MR(BC1.SMR.25.slopes,density=1000)

##Calculate Maximum Metabolic Rates from SMR Data
BC1.SMR.MMR.21.slopes<-extract.slope(BC1.SMR.21.clean,
                                 method="max",
                                 n.slope=1)

BC1.SMR.MMR.21.2.slopes<-extract.slope(BC1.SMR.21.2.clean,
                                 method="max",
                                 n.slope=1)

BC1.SMR.MMR.23.slopes<-extract.slope(BC1.SMR.23.clean,
                                 method="max",
                                 n.slope=1)

BC1.SMR.MMR.25.slopes<-extract.slope(BC1.SMR.25.clean,
                                 method="max",
                                 n.slope=1)

##Calculate Maximum Metabolic Rates
BC1.SMR.MMR.21<-calculate.MR(BC1.SMR.MMR.21.slopes,density=1000)
BC1.SMR.MMR.21.2<-calculate.MR(BC1.SMR.MMR.21.2.slopes,density=1000)
BC1.SMR.MMR.23<-calculate.MR(BC1.SMR.MMR.23.slopes,density=1000)
BC1.SMR.MMR.25<-calculate.MR(BC1.SMR.MMR.25.slopes,density=1000)


##Combine Data into one dataframe
BC1.SMR.Acute<-rbind.data.frame(BC1.SMR.21,BC1.SMR.21.2,BC1.SMR.23,BC1.SMR.25)
BC1.SMR.MMR.Acute<-rbind.data.frame(BC1.SMR.MMR.21,BC1.SMR.MMR.21.2,BC1.SMR.MMR.23,BC1.SMR.MMR.25)

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

BC1.SMR.MMR.Acute<-cbind(BC1.SMR.MMR.Acute,Trial,Temperature,Event)
BC1.SMR.Acute<-cbind(BC1.SMR.Acute,Temperature, Event)

########################################################################
##Save out CSVs
setwd("~/Desktop/Dissertation/Chapter 2")

write.csv(BC1.SMR.Acute,"LP BC Acute SMR EXP 1.csv")
write.csv(BC1.SMR.MMR.Acute,"LP BC Acute SMR.MMR EXP 1.csv")

#######################################################################
#Maximum Metabolic Rate
########################################################################

#Convert %Air Saturation to mg/L for MMR measurements and Empty chamber runs
#First Run of Acute MMR Upper Pennisula
convert.rMR('MMR_21_BC_8.4.2021_raw.txt',
            'Converted MMR 21 08.04.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_21.2_BC_8.7.2021_raw.txt',
            'Converted MMR 21.2 08.07.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_23_BC_8.5.2021_raw.txt',
            'Converted MMR 23 08.05.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_25_BC_8.6.2021_raw.txt',
            'Converted MMR 25 08.06.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw MMR Data
BC1.MMR.21.raw<-import.meas('Converted MMR 21 08.04.2021.txt', 
                            info.data=info_21.04, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
BC1.MMR.21.2.raw<-import.meas('Converted MMR 21.2 08.07.2021.txt', 
                              info.data=info_21.2.07, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
BC1.MMR.23.raw<-import.meas('Converted MMR 23 08.05.2021.txt', 
                            info.data=info_23.05, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
BC1.MMR.25.raw<-import.meas('Converted MMR 25 08.06.2021.txt', 
                            info.data=info_25.06, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)

#Clean Raw MMR data with background respiration measurments
BC1.MMR.21.clean<-correct.meas(info.data=info_21.04,
                               pre.data=BC1.pre,
                               post.data=BC1.post,
                               meas.data=BC1.MMR.21.raw,
                               method="exponential")
BC1.MMR.21.2.clean<-correct.meas(info.data=info_21.2.07,
                                 pre.data=BC1.pre,
                                 post.data=BC1.post, 
                                 meas.data=BC1.MMR.21.2.raw,
                                 method="exponential")
BC1.MMR.23.clean<-correct.meas(info.data=info_23.05,
                               pre.data=BC1.pre,
                               post.data=BC1.post, 
                               meas.data=BC1.MMR.23.raw,
                               method="exponential")
BC1.MMR.25.clean<-correct.meas(info.data=info_25.06,
                               pre.data=BC1.pre,
                               post.data=BC1.post,
                               meas.data=BC1.MMR.25.raw,
                               method="exponential")

#Extract slopes
BC1.MMR.21.slopes<-extract.slope(BC1.MMR.21.clean,
                                 method="max",
                                 n.slope=1)
BC1.MMR.21.2.slopes<-extract.slope(BC1.MMR.21.2.clean,
                                   method="max",
                                   n.slope=1)
BC1.MMR.23.slopes<-extract.slope(BC1.MMR.23.clean,
                                 method="max",
                                 n.slope=1)
BC1.MMR.25.slopes<-extract.slope(BC1.MMR.25.clean,
                                 method="max",
                                 n.slope=1)

##Calculate Metabolic Rates
BC1.MMR.21<-calculate.MR(BC1.MMR.21.slopes,density=1000)
BC1.MMR.21.2<-calculate.MR(BC1.MMR.21.2.slopes,density=1000)
BC1.MMR.23<-calculate.MR(BC1.MMR.23.slopes,density=1000)
BC1.MMR.25<-calculate.MR(BC1.MMR.25.slopes,density=1000)


##Combine into one Datafram
BC1.MMR.Acute<-rbind.data.frame(BC1.MMR.21,BC1.MMR.21.2,BC1.MMR.23,BC1.MMR.25)

##Add columns
Trial<-c("MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR","MMR","MMR")
BC1.MMR.Acute<-cbind(BC1.MMR.Acute,Trial,Temperature,Event)

#################3
##Save out as CSV

setwd("~/Desktop/Dissertation/Chapter 2")

write.csv(BC1.MMR.Acute,"LP BC Acute MMR EXP 1.csv")

#######################################################################
#Data Analysis and grapahs
#######################################################################

BC1.MMR<-rbind(BC1.SMR.MMR.Acute,BC1.MMR.Acute)
BC1.MMR<-read.csv("LP BC MMR BEST EXP 1.csv")
BC1.MMR$Temperature<-factor(Temperature)
BC1.SMR.MMR.Acute$Temperature<-factor(Temperature)
BC1.MMR.Acute$Temperature<-factor(Temperature)
MMRC$Temperature<-factor(Temperature)
MMRC<-as.data.frame(MMRC)
m<-ggplot() 
m<-m+geom_point(data=BC1.SMR.MMR.Acute,aes(Temperature, MR.mass/60, group=Ind,color=Ind)) + ggtitle("BC1 MMR, Comparison")+theme_classic()
m<-m+geom_point(data=BC1.MMR.Acute,aes(Temperature, MR.mass/60, group=Ind,color=Ind))
m<-m+ylab("MO2 (mgO2/(kg min))")
m

s<-ggplot(BC2.SMR,aes(Trial, MR.mass/60, group=Ind,color=Ind)) 
s<-s+geom_line() + geom_point() + ggtitle("BC2 SMR, Comparison")+theme_classic()
s<-s+ylab("MO2 (mgO2/(kg min))")


grid.arrange(
        ggplot(SMR.Acute,aes(x=Temperature,y=MR.mass/60),group=Stock,color=Stock) + 
                geom_boxplot() +geom_boxplot(data=BC1.MMR, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
                ylab("MO2 (mgO2/(kg min))") + ggtitle("BC SMR,Best MMR ") + theme(legend.position="none"),
        ggplot(BC1.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
                geom_boxplot() +geom_boxplot(data=BC1.MMR.Acute, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
                ylab("MO2 (mgO2/(kg min))") + ggtitle("BC SMR, MMR ") + theme(legend.position="none"),
        ggplot(BC1.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
                geom_boxplot() +geom_boxplot(data=BC1.SMR.MMR.Acute, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
                ylab("MO2 (mgO2/(kg min))") + ggtitle("BC SMR, SMR-MMR ") + theme(legend.position="none"),
        
        ncol=3
)

BC1.SMR.Acute$Stock <-  "LP"
LF3.SMR.Acute$Stock<-"UP"

SMR.Acute<-rbind(BC1.SMR.Acute,LF3.SMR.Acute)

BC1.MMR$Stock <- "LP"
LF3.MMR$Stock<-"UP"
BC1.MMR.Acute$Stock<- "LP"
LF3.MMR.Acute$Stock<-"UP"
MMR.Acute<-rbind(BC1.MMR,LF3.MMR)
BC1.SMR.MMR.Acute$Stock<-"LP"
MMR<-rbind(BC1.MMR.Acute,BC1.SMR.MMR.Acute)
spring_data$Temp.factor<-factor(spring_data$Temp.factor)
s<-ggplot(spring_data,aes(x=Temp.factor,y=SMR/60,color=Stock))+geom_boxplot()+geom_boxplot(data=spring_data, aes(x=Temp.factor,y=MMR/60,color=Stock))
s


grid.arrange(
        ggplot(spring_data,aes(x=as.factor(Temperature),y=SMR/60,group=Stock,color=Stock)) + geom_boxplot()
        +geom_boxplot(data=spring_data, aes(x=as.factor(Temperature),y=MMR/60,color=Stock))
        +ylab("MO2 (mgO2/(kg min))")+ylim(3,14) +xlab("Temperature")+ ggtitle("BC vs LF SMR,Global MMR ") + theme(legend.position="none"),
        ggplot(SMR.Acute,aes(x=as.factor(Temperature),y=MR.mass/60,color=Stock)) + 
                geom_boxplot() +geom_boxplot(data=MMR, aes(x=as.factor(Temperature),y=MR.mass/60,color=Stock))+ 
                ylab("MO2 (mgO2/(kg min))")+ylim(3,14) +xlab("Temperature") + ggtitle("BC vs LF SMR, MMR ") + theme(),
        ncol=2
)

grid.arrange(
        ggplot(subset(BC1.MMR,Event%in%c(1)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, BC 21") + theme(legend.position="none"),
        ggplot(subset(BC1.MMR,Event%in%c(4)), 
               aes(x = Trial, y = MR.mass,  group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, BC 21.2") + theme(legend.position="none"),
        ggplot(subset(BC1.MMR,Event%in%c(2)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, BC 23") + theme(legend.position="none"),
        ggplot(subset(BC1.MMR,Event%in%c(3)), 
               aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, BC 25") + theme(legend.position="none"),
        
        ggplot(subset(BC2.MMR,Event%in%c(1)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, BC 25") + theme(legend.position="none"),
        ggplot(subset(BC2.MMR,Event%in%c(4)), 
               aes(x = Trial, y = MR.mass,  group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, BC 25.2") + theme(legend.position="none"),
        ggplot(subset(BC2.MMR,Event%in%c(2)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, BC 23") + theme(legend.position="none"),
        ggplot(subset(BC2.MMR,Event%in%c(3)), 
               aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, BC 21") + theme(legend.position="none"),
        
        ggplot(subset(LF3.MMR,Event%in%c(1)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF.1 21") + theme(legend.position="none"),
        ggplot(subset(LF3.MMR,Event%in%c(4)), 
               aes(x = Trial, y = MR.mass,  group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF.1 21.2") + theme(legend.position="none"),
        ggplot(subset(LF3.MMR,Event%in%c(2)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF.1 23") + theme(legend.position="none"),
        ggplot(subset(LF3.MMR,Event%in%c(3)), 
               aes(x = Trial, y = MR.mass,  group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF.1 25") + theme(legend.position="none"),
        
        ggplot(subset(LF4.MMR,Event%in%c(1)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF.2 25") + theme(legend.position="none"),
        ggplot(subset(LF4.MMR,Event%in%c(4)), 
               aes(x = Trial, y = MR.mass,  group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF.2 25.2") + theme(legend.position="none"),
        ggplot(subset(LF4.MMR,Event%in%c(2)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF.2 23") + theme(legend.position="none"),
        ggplot(subset(LF4.MMR,Event%in%c(3)), 
               aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF.2 21") + theme(legend.position="none"),
        ncol=4
)





