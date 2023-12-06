library(FishResp)
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_7_26_2021")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used
info_21<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(3.53,2.85,1.78,2.56,2.6,2.91,3.19,2.11), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_23<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(3.53,2.85,0,2.56,0,2.91,3.19,2.11), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_25<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                                Mass=c(0,2.85,0,2.56,0,2.91,3.19,2.11), 
                                Volume = c(78,78,78,78,78,78,78,78), 
                                DO.unit = "mg/L")
info_21.2<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                        Mass=c(0,2.85,0,2.56,0,2.91,3.19,2.11), 
                        Volume = c(78,78,78,78,78,78,78,78), 
                        DO.unit = "mg/L")

#Convert blank runs used for background respiration
convert.rMR('empty_pre2_7_26_2021_raw.txt', 
            'Converted Empty_Pre 07.26.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_7.30.2021_raw.txt',
            'Converted Empty_Post 07.30.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")


#Import Background Respiration Files
LF1.pre<-import.test('Converted Empty_Pre 07.26.2021.txt', 
                     info.data=info_08, 
                     logger= "AutoResp",
                     n.chamber=8,
                     plot.temperature= F,
                     plot.oxygen = F)
LF1.post<-import.test('Converted Empty_Post 07.30.2021.txt',
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
convert.rMR('SMR_LF_21_7_26_2021_raw.txt',
            'Converted SMR 21 07.26.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_21.2_LF_7.29.2021_raw.txt',
            'Converted SMR 21.2 07.29.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_23_LF_7_27_2021_raw.txt',
            'Converted SMR 23 07.27.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_25_LF_7.28.2021_raw.txt',
            'Converted SMR 25 07.28.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw SMR Data
LF1.SMR.21.raw<-import.meas('Converted SMR 21 07.26.2021.txt', 
                        info.data=info_21, 
                        logger="AutoResp", 
                        n.chamber=8,
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)
LF1.SMR.21.2.raw<-import.meas('Converted SMR 21.2 07.29.2021.txt', 
                        info.data=info_21.2, 
                        n.chamber=8,
                        logger="AutoResp", 
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)
LF1.SMR.23.raw<-import.meas('Converted SMR 23 07.27.2021.txt', 
                          info.data=info_23, 
                          logger="AutoResp", 
                          n.chamber=8,
                          date.format="MDY",
                          plot.temperature = F,
                          plot.oxygen = F)
LF1.SMR.25.raw<-import.meas('Converted SMR 25 07.28.2021.txt', 
                        info.data=info_25, 
                        logger="AutoResp", 
                        n.chamber=8,
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)

#Clean Raw SMR data with background respiration measurments
LF1.SMR.21.clean<-correct.meas(info.data=info_21,
                               pre.data=LF1.pre,
                               post.data=LF1.post,
                           meas.data=LF1.SMR.21.raw,
                           method="exponential")
LF1.SMR.21.2.clean<-correct.meas(info.data=info_21.2,
                           pre.data=LF1.pre,
                           post.data=LF1.post, 
                           meas.data=LF1.SMR.21.2.raw,
                           method="exponential")
LF1.SMR.23.clean<-correct.meas(info.data=info_23,
                               pre.data=LF1.pre,
                               post.data=LF1.post, 
                             meas.data=LF1.SMR.23.raw,
                             method="exponential")
LF1.SMR.25.clean<-correct.meas(info.data=info_25,
                               pre.data=LF1.pre,
                               post.data=LF1.post,
                           meas.data=LF1.SMR.25.raw,
                           method="exponential")

#Extract slopes
LF1.SMR.21.slopes<-extract.slope(LF1.SMR.21.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.25)
LF1.SMR.21.2.slopes<-extract.slope(LF1.SMR.21.2.clean,
                             method="calcSMR.quant",
                             r2=0.5,
                             p=0.25)
LF1.SMR.23.slopes<-extract.slope(LF1.SMR.23.clean,
                               method="calcSMR.quant",
                               r2=0.5,
                               p=0.25)
LF1.SMR.25.slopes<-extract.slope(LF1.SMR.25.clean,
                             method="calcSMR.quant",
                             r2=0.55,
                             p=0.25)

##Calculate Metabolic Rates
LF1.SMR.21<-calculate.MR(LF1.SMR.21.slopes,density=1000)
LF1.SMR.21.2<-calculate.MR(LF1.SMR.21.2.slopes,density=1000)
LF1.SMR.23<-calculate.MR(LF1.SMR.23.slopes,density=1000)
LF1.SMR.25<-calculate.MR(LF1.SMR.25.slopes,density=1000)

#Combine Data into one dataframe
LF1.Acute<-rbind.data.frame(LF1.SMR.21,LF1.SMR.21.2,LF1.SMR.23,LF1.SMR.25)
Temperature<-c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,23,23,23,23,23,23,23,23,25,25,25,25,25,25,25,25)
LF1.Accute<-cbind(LF1.Acute, Temperature)
setwd("~/Desktop/Dissertation/Chapter 2")
write.csv(LF1.Accute,"UP LF Acute SMR EXP 1")

#######################################################################
#Maximum Metabolic Rate
########################################################################

#Convert %Air Saturation to mg/L for MMR measurements and Empty chamber runs
#First Run of Acute MMR Upper Pennisula
convert.rMR('MMR_UP_21_7_27_2021_raw.txt',
            'Converted MMR 21 07.27.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_21.2_LF_7.30.2021_raw.txt',
            'Converted MMR 21.2 07.30.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_23_LF_7.28.2021_raw.txt',
            'Converted MMR 23 07.28.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_25_LF_7.29.2021_raw.txt',
            'Converted MMR 25 07.29.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw MMR Data
LF1.MMR.21.raw<-import.meas('Converted MMR 21 07.27.2021.txt', 
                            info.data=info_21, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
LF1.MMR.21.2.raw<-import.meas('Converted MMR 21.2 07.30.2021.txt', 
                              info.data=info_21.2, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)
LF1.MMR.23.raw<-import.meas('Converted MMR 23 07.28.2021.txt', 
                            info.data=info_23, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
LF1.MMR.25.raw<-import.meas('Converted MMR 25 07.29.2021.txt', 
                            info.data=info_25, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)

#Clean Raw MMR data with background respiration measurments
LF1.MMR.21.clean<-correct.meas(info.data=info_21,
                               pre.data=LF1.pre,
                               post.data=LF1.post,
                               meas.data=LF1.MMR.21.raw,
                               method="exponential")
LF1.MMR.21.2.clean<-correct.meas(info.data=info_21.2,
                                 pre.data=LF1.pre,
                                 post.data=LF1.post, 
                                 meas.data=LF1.MMR.21.2.raw,
                                 method="exponential")
LF1.MMR.23.clean<-correct.meas(info.data=info_23,
                               pre.data=LF1.pre,
                               post.data=LF1.post, 
                               meas.data=LF1.MMR.23.raw,
                               method="exponential")
LF1.MMR.25.clean<-correct.meas(info.data=info_25,
                               pre.data=LF1.pre,
                               post.data=LF1.post,
                               meas.data=LF1.MMR.25.raw,
                               method="exponential")

#Extract slopes
LF1.MMR.21.slopes<-extract.slope(LF1.MMR.21.clean,
                                 method="max",
                                 n.slope=1)
LF1.MMR.21.2.slopes<-extract.slope(LF1.MMR.21.2.clean,
                                   method="max",
                                   n.slope=1)
LF1.MMR.23.slopes<-extract.slope(LF1.MMR.23.clean,
                                 method="max",
                                 n.slope=1)
LF1.MMR.25.slopes<-extract.slope(LF1.MMR.25.clean,
                                 method="max",
                                 n.slope=1)

##Calculate Metabolic Rates
LF1.MMR.21<-calculate.MR(LF1.MMR.21.slopes,density=1000)
LF1.MMR.21.2<-calculate.MR(LF1.MMR.21.2.slopes,density=1000)
LF1.MMR.23<-calculate.MR(LF1.MMR.23.slopes,density=1000)
LF1.MMR.25<-calculate.MR(LF1.MMR.25.slopes,density=1000)








#######################################################################
#Data Analysis and grapahs
#######################################################################


plot((UP1.Accute$MR.mass/1000)~UP1.Accute$Temperature,
     yaxt="none",
     ylab="Mass specific Metabolic Rate g/L/h",
     xlab= "Temperaturec(C)",
     main="Standard Metabolic Rates", 
     ylim=c(0,.5),
     xlim=c(21,25),
     cex=.75, col=1, pch=19)
axis(2,seq(0,2,.1))





