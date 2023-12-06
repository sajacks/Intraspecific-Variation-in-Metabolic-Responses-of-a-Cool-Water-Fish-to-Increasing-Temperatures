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
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_10.31.2021")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used
info_21.31<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                        Mass=c(19.38,13.78,15.84,16.23,15.06,15.98,11.55,13.24), 
                        Volume = c(317,317,317,317,317,317,317,317), 
                        DO.unit = "mg/L")
info_21.01<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(18.9,13.74,15.56,15.81,15.92,11.32,13.05), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_23.02<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(14.55,10.82,12.97,7.73,11.53,15.58,12.78,8.83),
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_23.03<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(14.66,10.81,12.72,7.49,11.38,15.48,13.69,8.53), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_25.04<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(17.61,18.05,13.32,21.4,12.42,15.72,12.03,12.36), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_25.05<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(17.25,17.8,13,21,12.24,15.31,11.86,12.19), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")

#Convert blank runs used for background respiration
convert.rMR('empty_post_pre_10.30.2021_raw.txt', 
            'Converted Empty_Pre 10.30.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_11.1.2021_raw.txt',
            'Converted Empty_Post 11.01.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_pre_11.1.2021_raw.txt',
            'Converted Empty_pre 11.01.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_pre_11.3.2021_raw.txt',
            'Converted Empty_Post_Pre 11.03.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('empty_post_11.5.2021_raw.txt',
            'Converted Empty_Post 11.05.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")



#Import Background Respiration Files
DB21.pre<-import.test('Converted Empty_Pre 10.30.2021.txt', 
                      info.data=info_21.31, 
                      logger= "AutoResp",
                      n.chamber=8,
                      plot.temperature= F,
                      plot.oxygen = F)
DB21.post<-import.test('Converted Empty_Post 11.01.2021.txt',
                       info.data=info_21.01, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)
DB23.pre<-import.test('Converted Empty_Pre 11.01.2021.txt',
                      info.data=info_23.02, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
DB23.post<-import.test('Converted Empty_Post_Pre 11.03.2021.txt',
                       info.data=info_23.03, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)

DB25.pre<-import.test('Converted Empty_Post_Pre 11.03.2021.txt',
                      info.data=info_25.04, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
DB25.post<-import.test('Converted Empty_post 11.05.2021.txt',
                       info.data=info_25.05, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)


#######################################################################
#Standard Metabolic Rate
#######################################################################
#Convert %Air Saturation to mg/L for SMR measurements and Empty chamber runs
#First Run of Acute SMR Lower Pennisula from 2019 Fall fingerling data
convert.rMR('SMR_21_DB_10.31.2021_raw.txt',
            'Converted SMR 21 10.31.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_23_DB_11.2.2021_raw.txt',
            'Converted SMR 23 11.02.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_25_DB_11.4.2021_raw.txt',
            'Converted SMR 25 11.04.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw SMR Data
DB.SMR.21.raw<-import.meas('Converted SMR 21 10.31.2021.txt', 
                              info.data=info_21.31, 
                              logger="AutoResp", 
                              n.chamber=8,
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
DB.SMR.23.raw<-import.meas('Converted SMR 23 11.02.2021.txt', 
                              info.data=info_23.02, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
DB.SMR.25.raw<-import.meas('Converted SMR 25 11.04.2021.txt', 
                              info.data=info_25.04, 
                              logger="AutoResp", 
                              n.chamber=8,
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)


#Clean Raw SMR data with background respiration measurments
DB.SMR.21.clean<-correct.meas(info.data=info_21.31,
                              pre.data=DB21.pre,
                              post.data=DB21.post,
                              meas.data=DB.SMR.21.raw,
                              method="exponential")
DB.SMR.23.clean<-correct.meas(info.data=info_23.02,
                              pre.data=DB23.pre,
                              post.data=DB23.post, 
                              meas.data=DB.SMR.23.raw,
                              method="exponential")
DB.SMR.25.clean<-correct.meas(info.data=info_25.04,
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


write.csv(DB.SMR.Acclimated,"DB Acclimated SMR EXP 2.csv")
write.csv(DB.SMR.MMR.Acclimated,"DB Acclimated SMR-MMR EXP 2.csv")


#######################################################################
#Maximum Metabolic Rate
########################################################################

#Convert %Air Saturation to mg/L for MMR measurements and Empty chamber runs
#First Run of Acute MMR Upper Pennisula
convert.rMR('MMR_21_DB_11.1.2021_raw.txt',
            'Converted MMR 21 11.01.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_23_DB_11.3.2021_raw.txt',
            'Converted MMR 23 11.03.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_25_DB_11.5.2021_raw.txt',
            'Converted MMR 25 11.05.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw MMR Data
DB.MMR.21.raw<-import.meas('Converted MMR 21 11.01.2021.txt', 
                            info.data=info_21.01, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
DB.MMR.23.raw<-import.meas('Converted MMR 23 11.03.2021.txt', 
                            info.data=info_23.03, 
                            n.chamber=8,
                            logger="AutoResp", 
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
DB.MMR.25.raw<-import.meas('Converted MMR 25 11.05.2021.txt', 
                            info.data=info_25.05, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)

#Clean Raw MMR data with background respiration measurments
DB.MMR.21.clean<-correct.meas(info.data=info_21.01,
                              pre.data=DB21.pre,
                              post.data=DB21.post,
                              meas.data=DB.MMR.21.raw,
                              method="exponential")
DB.MMR.23.clean<-correct.meas(info.data=info_23.03,
                              pre.data=DB23.pre,
                              post.data=DB23.post, 
                              meas.data=DB.MMR.23.raw,
                              method="exponential")
DB.MMR.25.clean<-correct.meas(info.data=info_25.05,
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


write.csv(DB.MMR.Acclimated,"DB Acclimated MMR EXP 2.csv")

##Determine the highest MMR between Chase and MMR and Save as BEST file
DB.MMR.BEST<-NULL
for(i in 1:length(DB.SMR.MMR.Acclimated$MR.mass)){
  if(DB.SMR.MMR.Acclimated$MR.mass[i]>DB.MMR.Acclimated$MR.mass[i]){
    DB.MMR.BEST<-rbind(DB.MMR.BEST,DB.SMR.MMR.Acclimated[i,])
  }else {
    DB.MMR.BEST<-rbind(DB.MMR.BEST,DB.MMR.Acclimated[i,])
  }
}
write.csv(DB.MMR.BEST, "DB Acclimated MMR BEST EXP 2.csv")

#######################################################################
#Data Analysis and grapahs
#######################################################################

DB.MMR<-rbind(DB.SMR.MMR.Acclimated,DB.MMR.Acclimated)

grid.arrange(
  ggplot(subset(DB.MMR,Event%in%c(1)), 
         aes(x = Trial, y = MR.mass, group=Ind,color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, DB 21") + theme(legend.position="none")+ylim(200,800),
  ggplot(subset(DB.MMR,Event%in%c(2)), 
         aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, DB 23") + theme(legend.position="none")+ylim(200,800),
  ggplot(subset(DB.MMR,Event%in%c(3)), 
         aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, DB 25") + theme(legend.position="none")+ylim(200,800),
  
  ncol=3
)

