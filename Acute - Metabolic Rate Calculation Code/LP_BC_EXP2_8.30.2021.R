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
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_8.30.2021")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used
info_25.31<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(3.88,4.19,4.56,4.33,3.99,4.2,4.79,3.99), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_25.01<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(3.77,4.07,4.46,4.20,3.82,4.08,4.73,3.85), 
                       Volume = c(78,78,78,78,78,78,78,78), 
                       DO.unit = "mg/L")
info_23.02<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(3.69,4.14,4.29,4.14,3.81,3.96,4.64,3.66), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_21.03<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                                Mass=c(3.61,4.02,4.26,4.02,3.76,3.9,4.64,3.7), 
                                Volume = c(78,78,78,78,78,78,78,78), 
                                DO.unit = "mg/L")
info_25.2.04<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                        Mass=c(3.56,3.86,4.23,3.90,3.7,3.79,4.47,3.59), 
                        Volume = c(78,78,78,78,78,78,78,78), 
                        DO.unit = "mg/L")

#Convert blank runs used for background respiration
convert.rMR('empty_pre_8.31.2021_raw.txt', 
            'Converted Empty_Pre 08.31.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_9.4.2021_raw.txt',
            'Converted Empty_Post 09.04.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")


#Import Background Respiration Files
BC2.pre<-import.test('Converted Empty_Pre 08.31.2021.txt', 
                     info.data=info_21.24, 
                     logger= "AutoResp",
                     n.chamber=8,
                     plot.temperature= F,
                     plot.oxygen = F)
BC2.post<-import.test('Converted Empty_Post 09.04.2021.txt',
                      info.data=info_21.24, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)

#######################################################################
#Standard Metabolic Rate
#######################################################################
#Convert %Air Saturation to mg/L for SMR measurements and Empty chamber runs
#First Run of Accute SMR Upper Pennisula
convert.rMR('SMR_25_BC_8.31.2021_raw.txt',
            'Converted SMR 25 08.31.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_23_BC_9.1.2021_raw.txt',
            'Converted SMR 23 09.01.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_21_BC_9.2.2021_raw.txt',
            'Converted SMR 21 09.02.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_25.2_BC_9.3.2021_raw.txt',
            'Converted SMR 25.2 09.03.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
#Import Raw SMR Data
BC2.SMR.25.raw<-import.meas('Converted SMR 25 08.31.2021.txt', 
                        info.data=info_25.31, 
                        logger="AutoResp", 
                        n.chamber=8,
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)
BC2.SMR.23.raw<-import.meas('Converted SMR 23 09.01.2021.txt', 
                          info.data=info_25.01, 
                          logger="AutoResp", 
                          n.chamber=8,
                          date.format="MDY",
                          plot.temperature = F,
                          plot.oxygen = F)
BC2.SMR.21.raw<-import.meas('Converted SMR 21 09.02.2021.txt',
                            info.data=info_23.02,
                            logger="AutoResp",
                            n.chamber=8,date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
BC2.SMR.25.2.raw<-import.meas('Converted SMR 25.2 09.03.2021.txt', 
                              info.data=info_21.03, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)

#Clean Raw SMR data with background respiration measurments
BC2.SMR.25.clean<-correct.meas(info.data=info_25.31,
                               pre.data=BC2.pre,
                               post.data=BC2.post,
                           meas.data=BC2.SMR.25.raw,
                           method="exponential")
BC2.SMR.23.clean<-correct.meas(info.data=info_25.01,
                               pre.data=BC2.pre,
                               post.data=BC2.post, 
                             meas.data=BC2.SMR.23.raw,
                             method="exponential")
BC2.SMR.21.clean<-correct.meas(info.data=info_23.02,
                               pre.data=BC2.pre,
                               post.data=BC2.post,
                           meas.data=BC2.SMR.21.raw,
                           method="exponential")
BC2.SMR.25.2.clean<-correct.meas(info.data=info_21.03,
                                 pre.data=BC2.pre,
                                 post.data=BC2.post, 
                                 meas.data=BC2.SMR.25.2.raw,
                                 method="exponential")

#Extract slopes for Standard Metabolic Rate
BC2.SMR.25.slopes<-extract.slope(BC2.SMR.25.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)

BC2.SMR.23.slopes<-extract.slope(BC2.SMR.23.clean,
                               method="calcSMR.quant",
                               r2=0.95,
                               p=0.2)
BC2.SMR.21.slopes<-extract.slope(BC2.SMR.21.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)
BC2.SMR.25.2.slopes<-extract.slope(BC2.SMR.25.2.clean,
                                   method="calcSMR.quant",
                                   r2=0.95,
                                   p=0.2)

##Extract slopes for Maximum Metabolic Rate during SMR Measurements
BC2.SMR.MMR.25.slopes<-extract.slope(BC2.SMR.25.clean,
                                 method="max",
                                 n.slope=1)
BC2.SMR.MMR.23.slopes<-extract.slope(BC2.SMR.23.clean,
                                 method="max",
                                 n.slope=1)
BC2.SMR.MMR.21.slopes<-extract.slope(BC2.SMR.21.clean,
                                 method="max",
                                 n.slope=1)
BC2.SMR.MMR.25.2.slopes<-extract.slope(BC2.SMR.25.2.clean,
                                       method="max",
                                       n.slope=1)



##Calculate Standard Metabolic Rates
BC2.SMR.25<-calculate.MR(BC2.SMR.25.slopes,density=1000)
BC2.SMR.23<-calculate.MR(BC2.SMR.23.slopes,density=1000)
BC2.SMR.21<-calculate.MR(BC2.SMR.21.slopes,density=1000)
BC2.SMR.25.2<-calculate.MR(BC2.SMR.25.2.slopes,density=1000)

##Calculate Maximum Metabolic Rates 
BC2.SMR.MMR.25<-calculate.MR(BC2.SMR.MMR.25.slopes,density=1000)
BC2.SMR.MMR.23<-calculate.MR(BC2.SMR.MMR.23.slopes,density=1000)
BC2.SMR.MMR.21<-calculate.MR(BC2.SMR.MMR.21.slopes,density=1000)
BC2.SMR.MMR.25.2<-calculate.MR(BC2.SMR.MMR.25.2.slopes,density=1000)


##Combine Data into one Dataframe
BC2.SMR.Acute<-rbind.data.frame(BC2.SMR.25,BC2.SMR.23,BC2.SMR.21,BC2.SMR.25.2)
BC2.SMR.MMR.Acute<-rbind.data.frame(BC2.SMR.MMR.25,BC2.SMR.MMR.23,BC2.SMR.MMR.21,BC2.SMR.MMR.25.2)

##Add additional Columns
Temperature<-c(25,25,25,25,25,25,25,25,23,23,23,23,23,23,23,23,21,21,21,21,21,21,21,21,25,25,25,25,25,25,25,25)
Trial<-c("SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR","SMR","SMR")
Event<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4)

BC2.SMR.Acute<-cbind(BC2.SMR.Acute,Trial,Temperature,Event)
BC2.SMR.MMR.Acute<-cbind(BC2.SMR.MMR.Acute,Trial,Temperature,Event)

##############33
##Save out CSVs

write.csv(BC2.SMR.Acute,"LP BC Acute SMR EXP 2.csv")
write.csv(BC2.SMR.MMR.Acute,"LP BC Acute SMR-MMR EXP 2.csv")


#######################################################################
#Maximum Metabolic Rate
########################################################################

#Convert %Air Saturation to mg/L for MMR measurements and Empty chamber runs
#First Run of Acute MMR Upper Pennisula
convert.rMR('MMR_25_BC_9.1.2021_raw.txt',
            'Converted MMR 25 09.01.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_23_BC_9.2.2021_raw.txt',
            'Converted MMR 23 09.02.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_21_BC_9.3.2021_raw.txt',
            'Converted MMR 21 09.03.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_25_BC_9.4.2021_raw.txt',
            'Converted MMR 25.2 09.04.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw MMR Data
BC2.MMR.25.raw<-import.meas('Converted MMR 25 09.01.2021.txt', 
                            info.data=info_25.01, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
BC2.MMR.23.raw<-import.meas('Converted MMR 23 09.02.2021.txt', 
                            info.data=info_23.02, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
BC2.MMR.21.raw<-import.meas('Converted MMR 21 09.03.2021.txt', 
                            info.data=info_21.03, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
BC2.MMR.25.2.raw<-import.meas('Converted MMR 25.2 09.04.2021.txt', 
                              info.data=info_25.2.04, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)

#Clean Raw MMR data with background respiration measurments
BC2.MMR.25.clean<-correct.meas(info.data=info_25.01,
                               pre.data=BC2.pre,
                               post.data=BC2.post,
                               meas.data=BC2.MMR.25.raw,
                               method="exponential")
                                
BC2.MMR.23.clean<-correct.meas(info.data=info_23.02,
                               pre.data=BC2.pre,
                               post.data=BC2.post, 
                               meas.data=BC2.MMR.23.raw,
                               method="exponential")
BC2.MMR.21.clean<-correct.meas(info.data=info_21.03,
                               pre.data=BC2.pre,
                               post.data=BC2.post,
                               meas.data=BC2.MMR.21.raw,
                               method="exponential")
BC2.MMR.25.2.clean<-correct.meas(info.data=info_25.2.04,
                                 pre.data=BC2.pre,
                                 post.data=BC2.post, 
                                 meas.data=BC2.MMR.25.2.raw,
                                 method="exponential")

#Extract slopes
BC2.MMR.25.slopes<-extract.slope(BC2.MMR.25.clean,
                                 method="max",
                                 n.slope=1)
BC2.MMR.23.slopes<-extract.slope(BC2.MMR.23.clean,
                                 method="max",
                                 n.slope=1)
BC2.MMR.21.slopes<-extract.slope(BC2.MMR.21.clean,
                                 method="max",
                                 n.slope=1)
BC2.MMR.25.2.slopes<-extract.slope(BC2.MMR.25.2.clean,
                                   method="max",
                                   n.slope=1)

##Calculate Metabolic Rates
BC2.MMR.25<-calculate.MR(BC2.MMR.25.slopes,density=1000)
BC2.MMR.23<-calculate.MR(BC2.MMR.23.slopes,density=1000)
BC2.MMR.21<-calculate.MR(BC2.MMR.21.slopes,density=1000)
BC2.MMR.25.2<-calculate.MR(BC2.MMR.25.2.slopes,density=1000)


##Combine data into one Dataframe
BC2.MMR.Acute<-rbind.data.frame(BC2.MMR.25,BC2.MMR.23,BC2.MMR.21,BC2.MMR.25.2)

##Add additional columns
Trial<-c("MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR","MMR","MMR")
BC2.MMR.Acute<-cbind(BC2.MMR.Acute,Trial,Temperature,Event)

##Save data our as CSV


write.csv(BC2.MMR.Acute,"LP BC Acute MMR EXP 2.csv")

#######################################################################
#Data Analysis and grapahs
#######################################################################

BC2.MMR<-rbind(BC2.SMR.MMR.Acute,BC2.MMR.Acute)

BC2.MMR<-read.csv("LP BC Acute MMR BEST EXP 2.csv")
BC2.MMR$Temperature<-factor(Temperature)
BC2.SMR.MMR.Acute$Temperature<-factor(Temperature)
BC2.MMR.Acute$Temperature<-factor(Temperature)
BC2.SMR.Acute$Temperature<-factor(Temperature)




grid.arrange(
        ggplot(BC2.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
                geom_boxplot() +geom_boxplot(data=BC2.MMR, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
                ylab("MO2 (mgO2/(kg min))") + ggtitle("BC SMR,Best MMR ") + theme(legend.position="none"),
ggplot(BC2.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
        geom_boxplot() +geom_boxplot(data=BC2.MMR.Acute, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
        ylab("MO2 (mgO2/(kg min))") + ggtitle("BC SMR, MMR ") + theme(legend.position="none"),
ggplot(BC2.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
        geom_boxplot() +geom_boxplot(data=BC2.SMR.MMR.Acute, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
        ylab("MO2 (mgO2/(kg min))") + ggtitle("BC SMR, SMR-MMR ") + theme(legend.position="none"),
        
        ncol=3
)

grid.arrange(
        ggplot(subset(BC2.MMR,Event%in%c(1)), 
               aes(x = Trial, y = MR.mass, group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, BC 25") + theme(legend.position="none")+ylim(200,900),
        ggplot(subset(BC2.MMR,Event%in%c(4)), 
               aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, BC 25.2") + theme(legend.position="none")+ylim(200,900),
        ggplot(subset(BC2.MMR,Event%in%c(2)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, BC 23") + theme(legend.position="none")+ylim(200,900),
        ggplot(subset(BC2.MMR,Event%in%c(3)), 
               aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, BC 21") + theme(legend.position="none")+ylim(200,900),
        
        ncol=2
)
aMMR<-aov(MMRC$MR.mass~MMRC$Temperature*MMRC$Trial,data=MMRC)
Anova(aMMR)
TukeyHSD((aMMR))

as.factor(BC2.MMR$Temperature)


