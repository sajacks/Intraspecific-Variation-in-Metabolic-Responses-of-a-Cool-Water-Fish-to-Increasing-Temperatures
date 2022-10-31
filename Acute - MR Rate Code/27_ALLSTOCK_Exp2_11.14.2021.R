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
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_11.14.2021")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used
info_27.14<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                        Mass=c(18.16,16.05,21.85,13.35,13.9,20.63,18.83,13.067), 
                        Volume = c(317,317,317,317,317,317,317,317), 
                        DO.unit = "mg/L")
info_27.2.15<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(17.8,15.62,21.33,13.21,13.26,19.7,18.36,12.61), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_29.16<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(17.97,15.32,21,13.06,13.18,20.72,17.85,12.75), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_27.17<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(14.99,14.38,17.46,20.23,14.45,13.9,15.69,15.2), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_27.2.18<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                         Mass=c(14.51,13.95,17.16,19.28,14.12,13.55,15.14,14.82), 
                         Volume = c(317,317,317,317,317,317,317,317), 
                         DO.unit = "mg/L")
info_29.19<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(14.12,13.73,17.05,19.25,14.04,13.54,14.74,14.75), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_27.21<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(18.97,18.37,20.46,14.96,20.52,19.99,20.476,19.38), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")
info_27.2.22<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                         Mass=c(18.86,18.07,20.02,14.85,20.09,19.75,20.30,19.03), 
                         Volume = c(317,317,317,317,317,317,317,317), 
                         DO.unit = "mg/L")
info_29.23<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(18.3,17.49,20.07,14.43,19.83,19.64,19.75,18.81), 
                       Volume = c(317,317,317,317,317,317,317,317), 
                       DO.unit = "mg/L")

#Convert blank runs used for background respiration
convert.rMR('empty_pre2_11.13.2021_raw.txt', 
            'Converted Empty_Pre 11.13.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_pre_11.16.2021_raw.txt',
            'Converted Empty_post_pre 11.16.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_11.19.2021_raw.txt',
            'Converted Empty_Post 11.19.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")


convert.rMR('empty_pre_11.20.2021_raw.txt', 
            'Converted Empty_Pre 11.20.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_11.23.2021_raw.txt',
            'Converted Empty_post 11.23.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")



#Import Background Respiration Files
DB.pre<-import.test('Converted Empty_Pre 11.13.2021.txt', 
                      info.data=info_27.14, 
                      logger= "AutoResp",
                      n.chamber=8,
                      plot.temperature= F,
                      plot.oxygen = F)
DB.post<-import.test('Converted Empty_post_pre 11.16.2021.txt',
                       info.data=info_29.16, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)
BC.pre<-import.test('Converted Empty_post_pre 11.16.2021.txt',
                      info.data=info_27.17, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
BC.post<-import.test('Converted Empty_Post 11.19.2021.txt',
                       info.data=info_29.19, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)

LF.pre<-import.test('Converted Empty_Pre 11.20.2021.txt',
                      info.data=info_27.21, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)
LF.post<-import.test('Converted Empty_post 11.23.2021.txt',
                       info.data=info_29.23, 
                       logger="AutoResp",
                       n.chamber=8,
                       plot.temperature = F,
                       plot.oxygen = F)


#######################################################################
#Standard Metabolic Rate
#######################################################################
#Convert %Air Saturation to mg/L for SMR measurements and Empty chamber runs
#First Run of Acute SMR Lower Pennisula from 2019 Fall fingerling data
convert.rMR('SMR_27_DB_11.14.2021_raw.txt',
            'Converted SMR 27 11.14.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_29_DB_11.15.2021_raw.txt',
            'Converted SMR 29 11.15.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_27_BC_11.17.2021_raw.txt',
            'Converted SMR 27 11.17.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_29_BC_11.18.2021_raw.txt',
            'Converted SMR 29 11.18.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('SMR_27_LF_11.21.2021_raw.txt',
            'Converted SMR 27 11.21.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_29_LF_11.22.2021_raw.txt',
            'Converted SMR 29 11.22.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw SMR Data
DB.SMR.27.raw<-import.meas('Converted SMR 27 11.14.2021.txt', 
                           info.data=info_27.14, 
                           logger="AutoResp", 
                           n.chamber=8,
                           date.format="MDY",
                           plot.temperature = F,
                           plot.oxygen = F)
DB.SMR.29.raw<-import.meas('Converted SMR 29 11.15.2021.txt', 
                           info.data=info_27.2.15, 
                           n.chamber=8,
                           logger="AutoResp", 
                           date.format="MDY",
                           plot.temperature = F,
                           plot.oxygen = F)

BC.SMR.27.raw<-import.meas('Converted SMR 27 11.17.2021.txt', 
                              info.data=info_27.17, 
                              logger="AutoResp", 
                              n.chamber=8,
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
BC.SMR.29.raw<-import.meas('Converted SMR 29 11.18.2021.txt', 
                              info.data=info_27.2.18, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)

LF.SMR.27.raw<-import.meas('Converted SMR 27 11.21.2021.txt', 
                          info.data=info_27.21, 
                          logger="AutoResp", 
                          n.chamber=8,
                          date.format="MDY",
                          plot.temperature = F,
                          plot.oxygen = F)
LF.SMR.29.raw<-import.meas('Converted SMR 29 11.22.2021.txt', 
                           info.data=info_27.2.22, 
                           n.chamber=8,
                           logger="AutoResp", 
                           date.format="MDY",
                           plot.temperature = F,
                           plot.oxygen = F)


#Clean Raw SMR data with background respiration measurments
DB.SMR.27.clean<-correct.meas(info.data=info_27.14,
                              pre.data=DB.pre,
                              post.data=DB.post,
                              meas.data=DB.SMR.27.raw,
                              method="exponential")
DB.SMR.29.clean<-correct.meas(info.data=info_27.2.15,
                              pre.data=DB.pre,
                              post.data=DB.post, 
                              meas.data=DB.SMR.29.raw,
                              method="exponential")

BC.SMR.27.clean<-correct.meas(info.data=info_27.17,
                              pre.data=BC.pre,
                              post.data=BC.post,
                              meas.data=BC.SMR.27.raw,
                              method="exponential")
BC.SMR.29.clean<-correct.meas(info.data=info_27.2.18,
                              pre.data=BC.pre,
                              post.data=BC.post, 
                              meas.data=BC.SMR.29.raw,
                              method="exponential")

LF.SMR.27.clean<-correct.meas(info.data=info_27.21,
                              pre.data=LF.pre,
                              post.data=LF.post,
                              meas.data=LF.SMR.27.raw,
                              method="exponential")
LF.SMR.29.clean<-correct.meas(info.data=info_27.2.22,
                              pre.data=LF.pre,
                              post.data=LF.post, 
                              meas.data=LF.SMR.29.raw,
                              method="exponential")


#Extract slopes for Standard Metabolic Rate
DB.SMR.27.slopes<-extract.slope(DB.SMR.27.clean,
                                 method="calcSMR.quant",
                                 r2=0.9,
                                 p=0.2)
DB.SMR.29.slopes<-extract.slope(DB.SMR.29.clean,
                                 method="calcSMR.quant",
                                 r2=c(0.9),
                                 p=0.2)

BC.SMR.27.slopes<-extract.slope(BC.SMR.27.clean,
                                method="calcSMR.quant",
                                r2=0.9,
                                p=0.2)
BC.SMR.29.slopes<-extract.slope(BC.SMR.29.clean,
                                method="calcSMR.quant",
                                r2=c(0.9),
                                p=0.2)

LF.SMR.27.slopes<-extract.slope(LF.SMR.27.clean,
                                method="calcSMR.quant",
                                r2=0.9,
                                p=0.2)
LF.SMR.29.slopes<-extract.slope(LF.SMR.29.clean,
                                method="calcSMR.quant",
                                r2=c(0.9),
                                p=0.2)



##Extract slopes for Maximum Metabolic Rate during SMR Measurements
DB.SMR.MMR.27.slopes<-extract.slope(DB.SMR.27.clean,
                                     method="max",
                                     n.slope=1)
DB.SMR.MMR.29.slopes<-extract.slope(DB.SMR.29.clean,
                                     method="max",
                                     n.slope=1)

BC.SMR.MMR.27.slopes<-extract.slope(BC.SMR.27.clean,
                                    method="max",
                                    n.slope=1)
BC.SMR.MMR.29.slopes<-extract.slope(BC.SMR.29.clean,
                                    method="max",
                                    n.slope=1)

LF.SMR.MMR.27.slopes<-extract.slope(LF.SMR.27.clean,
                                    method="max",
                                    n.slope=1)
LF.SMR.MMR.29.slopes<-extract.slope(LF.SMR.29.clean,
                                    method="max",
                                    n.slope=1)




##Calculate Standard Metabolic Rates
DB.SMR.27<-calculate.MR(DB.SMR.27.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)
DB.SMR.29<-calculate.MR(DB.SMR.29.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)

BC.SMR.27<-calculate.MR(BC.SMR.27.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)
BC.SMR.29<-calculate.MR(BC.SMR.29.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)

LF.SMR.27<-calculate.MR(LF.SMR.27.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)
LF.SMR.29<-calculate.MR(LF.SMR.29.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)



##Calculate Maximum Metabolic Rates 
DB.SMR.MMR.27<-calculate.MR(DB.SMR.MMR.27.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)
DB.SMR.MMR.29<-calculate.MR(DB.SMR.MMR.29.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)

BC.SMR.MMR.27<-calculate.MR(BC.SMR.MMR.27.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)
BC.SMR.MMR.29<-calculate.MR(BC.SMR.MMR.29.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)

LF.SMR.MMR.27<-calculate.MR(LF.SMR.MMR.27.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)
LF.SMR.MMR.29<-calculate.MR(LF.SMR.MMR.29.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)



##Combine Data into one Dataframe
DB.SMR<-rbind.data.frame(DB.SMR.27,DB.SMR.29)
DB.SMR.MMR<-rbind.data.frame(DB.SMR.MMR.27,DB.SMR.MMR.29)

BC.SMR<-rbind.data.frame(BC.SMR.27,BC.SMR.29)
BC.SMR.MMR<-rbind.data.frame(BC.SMR.MMR.27,BC.SMR.MMR.29)

LF.SMR<-rbind.data.frame(LF.SMR.27,LF.SMR.29)
LF.SMR.MMR<-rbind.data.frame(LF.SMR.MMR.27,LF.SMR.MMR.29)

##Add additional Columns
Temperature<-c(27,27,27,27,27,27,27,27,29,29,29,29,29,29,29,29)
Trial<-c("SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR")
Event<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2)

DB.SMR<-cbind(DB.SMR,Trial,Temperature,Event)
DB.SMR.MMR<-cbind(DB.SMR.MMR,Trial,Temperature,Event)

BC.SMR<-cbind(BC.SMR,Trial,Temperature,Event)
BC.SMR.MMR<-cbind(BC.SMR.MMR,Trial,Temperature,Event)

LF.SMR<-cbind(LF.SMR,Trial,Temperature,Event)
LF.SMR.MMR<-cbind(LF.SMR.MMR,Trial,Temperature,Event)


##############
##Combine stocks into one dataframe for saving out
SMR.27<-rbind(DB.SMR,BC.SMR,LF.SMR)
SMR.MMR.27<-rbind(DB.SMR.MMR,BC.SMR.MMR,LF.SMR.MMR)

##############
##Save out CSVs


write.csv(SMR.27,"ALL_27_29 SMR Acute EXP 2.csv")
write.csv(SMR.MMR.27,"ALL_27_29 SMR.MMR Acute EXP 2.csv")


#######################################################################
#Maximum Metabolic Rate
########################################################################

#Convert %Air Saturation to mg/L for MMR measurements and Empty chamber runs
#First Run of Acute MMR Upper Pennisula
convert.rMR('MMR_27_DB_11.15.2021_raw.txt',
            'Converted MMR 27 11.15.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_29_DB_11.16.2021_raw.txt',
            'Converted MMR 29 11.16.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('MMR_27_BC_11.18.2021_raw.txt',
            'Converted MMR 27 11.18.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_29_BC_11.19.2021_raw.txt',
            'Converted MMR 29 11.19.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

convert.rMR('MMR_27_LF_11.22.2021_raw.txt',
            'Converted MMR 27 11.22.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_29_LF_11.23.2021_raw.txt',
            'Converted MMR 29 11.23.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")


#Import Raw MMR Data
DB.MMR.27.raw<-import.meas('Converted MMR 27 11.15.2021.txt', 
                            info.data=info_27.2.15, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
DB.MMR.29.raw<-import.meas('Converted MMR 29 11.16.2021.txt', 
                            info.data=info_29.16, 
                            n.chamber=8,
                            logger="AutoResp", 
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)

BC.MMR.27.raw<-import.meas('Converted MMR 27 11.18.2021.txt', 
                           info.data=info_27.2.18, 
                           logger="AutoResp", 
                           n.chamber=8,
                           date.format="MDY",
                           plot.temperature = F,
                           plot.oxygen = F)
BC.MMR.29.raw<-import.meas('Converted MMR 29 11.19.2021.txt', 
                           info.data=info_29.19, 
                           n.chamber=8,
                           logger="AutoResp", 
                           date.format="MDY",
                           plot.temperature = F,
                           plot.oxygen = F)

LF.MMR.27.raw<-import.meas('Converted MMR 27 11.22.2021.txt', 
                           info.data=info_27.2.22, 
                           logger="AutoResp", 
                           n.chamber=8,
                           date.format="MDY",
                           plot.temperature = F,
                           plot.oxygen = F)
LF.MMR.29.raw<-import.meas('Converted MMR 29 11.23.2021.txt', 
                           info.data=info_29.23, 
                           n.chamber=8,
                           logger="AutoResp", 
                           date.format="MDY",
                           plot.temperature = F,
                           plot.oxygen = F)


#Clean Raw MMR data with background respiration measurments
DB.MMR.27.clean<-correct.meas(info.data=info_27.2.15,
                              pre.data=DB.pre,
                              post.data=DB.post,
                              meas.data=DB.MMR.27.raw,
                              method="exponential")
DB.MMR.29.clean<-correct.meas(info.data=info_29.16,
                              pre.data=DB.pre,
                              post.data=DB.post, 
                              meas.data=DB.MMR.29.raw,
                              method="exponential")

BC.MMR.27.clean<-correct.meas(info.data=info_27.2.18,
                              pre.data=BC.pre,
                              post.data=BC.post,
                              meas.data=BC.MMR.27.raw,
                              method="exponential")
BC.MMR.29.clean<-correct.meas(info.data=info_29.19,
                              pre.data=BC.pre,
                              post.data=BC.post, 
                              meas.data=BC.MMR.29.raw,
                              method="exponential")

LF.MMR.27.clean<-correct.meas(info.data=info_27.2.22,
                              pre.data=LF.pre,
                              post.data=LF.post,
                              meas.data=LF.MMR.27.raw,
                              method="exponential")
LF.MMR.29.clean<-correct.meas(info.data=info_29.23,
                              pre.data=LF.pre,
                              post.data=LF.post, 
                              meas.data=LF.MMR.29.raw,
                              method="exponential")


#Extract slopes
DB.MMR.27.slopes<-extract.slope(DB.MMR.27.clean,
                                 method="max",
                                 r2=0.95,
                                 n.slope=1)
DB.MMR.29.slopes<-extract.slope(DB.MMR.29.clean,
                                 method="max",
                                 r2=0.95,
                                 n.slope=1)

BC.MMR.27.slopes<-extract.slope(BC.MMR.27.clean,
                                method="max",
                                r2=0.95,
                                n.slope=1)
BC.MMR.29.slopes<-extract.slope(BC.MMR.29.clean,
                                method="max",
                                r2=0.95,
                                n.slope=1)

LF.MMR.27.slopes<-extract.slope(LF.MMR.27.clean,
                                method="max",
                                r2=0.95,
                                n.slope=1)
LF.MMR.29.slopes<-extract.slope(LF.MMR.29.clean,
                                method="max",
                                r2=0.95,
                                n.slope=1)


##Calculate Metabolic Rates
DB.MMR.27<-calculate.MR(DB.MMR.27.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)
DB.MMR.29<-calculate.MR(DB.MMR.29.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)

BC.MMR.27<-calculate.MR(BC.MMR.27.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)
BC.MMR.29<-calculate.MR(BC.MMR.29.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)

LF.MMR.27<-calculate.MR(LF.MMR.27.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)
LF.MMR.29<-calculate.MR(LF.MMR.29.slopes,density=1000, plot.BR = F,plot.MR.abs = F,plot.MR.mass = F)



##Combine data into one Dataframe
DB.MMR<-rbind.data.frame(DB.MMR.27,DB.MMR.29)

BC.MMR<-rbind.data.frame(BC.MMR.27,BC.MMR.29)

LF.MMR<-rbind.data.frame(LF.MMR.27,LF.MMR.29)

##Add additional columns
Trial<-c("MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR")
Temperature<-c(27,27,27,27,27,27,27,27,29,29,29,29,29,29,29,29)
Event<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2)
DB.MMR<-cbind(DB.MMR,Trial,Temperature,Event)
BC.MMR<-cbind(BC.MMR,Trial,Temperature,Event)
LF.MMR<-cbind(LF.MMR,Trial,Temperature,Event)

##Combine into one dataframe
MMR.27<-rbind(DB.MMR,BC.MMR,LF.MMR)

##Save data our as CSV


write.csv(MMR.27,"ALL_27_29 Acute MMR EXP 2.csv")

##Determine the highest MMR between Chase and MMR and Save as BEST file
MMR.BEST<-NULL
length(MMR.27$MR.mass)
for(i in 1:length(MMR.27$MR.mass)){
  if(SMR.MMR.27$MR.mass[i]>MMR.27$MR.mass[i]){
    MMR.BEST<-rbind(MMR.BEST,SMR.MMR.27[i,])
  }else {
    MMR.BEST<-rbind(MMR.BEST,MMR.27[i,])
  }
}
write.csv(MMR.BEST, "ALL_27_29 Acute MMR BEST EXP 2.csv")

#######################################################################
#Data Analysis and grapahs
#######################################################################

MMR<-rbind(MMR.27,SMR.MMR.27)

grid.arrange(
  ggplot(subset(MMR,Event%in%c(1)), 
         aes(x = Trial, y = MR.mass, group=Ind,color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, 27") + theme(legend.position="none")+ylim(200,700),
  ggplot(subset(MMR,Event%in%c(2)), 
         aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
    geom_point()+geom_line() + ggtitle("MMR, 29") + theme(legend.position="none")+ylim(200,700),
  
  ncol=2
)





