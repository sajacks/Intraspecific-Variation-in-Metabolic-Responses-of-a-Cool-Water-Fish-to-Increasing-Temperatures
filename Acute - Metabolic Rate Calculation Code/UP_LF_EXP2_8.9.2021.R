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
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_8.9.2021")
#This experiment ended in the deaths of all eight fish. Post experiment background respiration was taken on after the fish were removed from the system. Likely an hour or so after death. 
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used
info_21<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(3.53,2.85,1.78,2.56,2.6,2.91,3.19,2.11), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_23<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(3.59,2.98,3.19,3.67,2.91,3.02,2.69,2.69), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_25<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                                Mass=c(), 
                                Volume = c(78,78,78,78,78,78,78,78), 
                                DO.unit = "mg/L")
info_21.2<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                        Mass=c(), 
                        Volume = c(78,78,78,78,78,78,78,78), 
                        DO.unit = "mg/L")

#Convert blank runs used for background respiration
convert.rMR('empty_pre_8.9.2021_raw.txt', 
            'Converted Empty_Pre 08.09.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_8.12.2021_raw.txt',
            'Converted Empty_Post 08.12.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")


#Import Background Respiration Files
LF2.pre<-import.test('Converted Empty_Pre 08.09.2021.txt', 
                     info.data=info_08, 
                     logger= "AutoResp",
                     n.chamber=8,
                     plot.temperature= F,
                     plot.oxygen = F)
LF2.post<-import.test('Converted Empty_Post 08.12.2021.txt',
                      info.data=info_08, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)

#######################################################################
#Standard Metabolic Rate
#######################################################################
#Convert %Air Saturation to mg/L for SMR measurements and Empty chamber runs
#First Run of Accute SMR Upper Pennisula
convert.rMR('SMR_23_LF_8.10.2021_raw.txt',
            'Converted SMR 23 08.10.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_23_LF_8.11.2021_raw.txt',
            'Converted MMR 23 08.11.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw SMR Data
LF2.SMR.23.raw<-import.meas('Converted SMR 23 08.10.2021.txt', 
                          info.data=info_23, 
                          logger="AutoResp", 
                          n.chamber=8,
                          date.format="MDY",
                          plot.temperature = F,
                          plot.oxygen = F)
LF2.MMR.23.raw<-import.meas('Converted MMR 23 08.11.2021.txt', 
                        info.data=info_25, 
                        logger="AutoResp", 
                        n.chamber=8,
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)

#Clean Raw SMR data with background respiration measurments

LF2.SMR.23.clean<-correct.meas(info.data=info_23,
                               pre.data=LF2.pre,
                               post.data=LF2.post, 
                             meas.data=LF2.SMR.23.raw,
                             method="exponential")
LF2.MMR.23.clean<-correct.meas(info.data=info_23,
                               pre.data=LF2.pre,
                               post.data=LF2.post, 
                               meas.data=LF2.MMR.23.raw,
                               method="exponential")


#Extract slopes

LF2.SMR.SMR.23.slopes<-extract.slope(LF2.SMR.23.clean,
                               method="calcSMR.quant",
                               r2=0.9,
                               p=0.2)
LF2.MMR.SMR.23.slopes<-extract.slope(LF2.MMR.23.clean,
                                 method="calcSMR.quant",
                                 r2=0.9,
                                 p=0.2)

LF2.SMR.MMR.23.slopes<-extract.slope(LF2.SMR.23.clean,
                                 method="max",
                                 n.slope=1)
LF2.MMR.MMR.23.slopes<-extract.slope(LF2.MMR.23.clean,
                                 method="max",
                                 n.slope=1)


##Calculate Metabolic Rates

LF2.SMR.SMR.23<-calculate.MR(LF2.SMR.SMR.23.slopes,density=1000)
LF2.MMR.SMR.23<-calculate.MR(LF2.MMR.SMR.23.slopes,density=1000)
LF2.SMR.MMR.23<-calculate.MR(LF2.SMR.MMR.23.slopes,density=1000)
LF2.MMR.MMR.23<-calculate.MR(LF2.MMR.MMR.23.slopes,density=1000)

#Combine Data into one dataframe
LF2.MMRC.Acute<-rbind.data.frame(LF2.SMR.MMR.23,LF2.MMR.MMR.23)
Trial<-c("SMR","SMR","SMR","SMR","SMR","SMR","SMR","SMR","MMR","MMR","MMR","MMR","MMR","MMR","MMR","MMR")
LF2.SMRC.Acute<-rbind.data.frame(LF2.SMR.SMR.23,LF2.MMR.SMR.23)

LF2.MMR<-cbind(LF2.MMRC.Acute,Trial)
LF2.SMR<-cbind(LF2.SMRC.Acute,Trial)


setwd("~/Desktop/")
write.csv(LF2.MMR,"UP_LF Acute MMR Comparison")
write.csv(LF2.SMR,"UP_LF Acute SMR Comparison")




#######################################################################
#Data Analysis and grapahs
#######################################################################


m<-ggplot(LF2.MMR,aes(Trial, MR.mass/60, group=Ind,color=Ind)) 
m<-m+geom_line() + geom_point() + ggtitle("LF2 MMR, Comparison")+theme_classic()
m<-m+ylab("MO2 (mgO2/(kg min))")
m

s<-ggplot(LF2.SMR,aes(Trial, MR.mass/60, group=Ind,color=Ind)) 
s<-s+geom_line() + geom_point() + ggtitle("LF2 SMR, Comparison")+theme_classic()
s<-s+ylab("MO2 (mgO2/(kg min))")
s


