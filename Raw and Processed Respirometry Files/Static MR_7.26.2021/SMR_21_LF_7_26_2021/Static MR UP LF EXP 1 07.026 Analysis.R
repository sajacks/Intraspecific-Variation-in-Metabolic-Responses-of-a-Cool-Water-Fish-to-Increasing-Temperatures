library(FishResp)
setwd("~/Desktop/MR data/Static MR_2019.12.08")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used
info_08<-input.info(ID=c(1,2,3), 
                    Mass=c(7.3,7.5,6.9), 
                    Volume = c(317,317,317), 
                    DO.unit = "mg/L")
#Convert %Air Saturation to mg/L for SMR measurements and Empty chamber runs
#First Run of Accute SMR Upper Pennisula
convert.rMR('SMR 21 2019.12.08_raw.txt',
            'Converted SMR 21 2019.12.08.txt',
            n.chamber=3,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR 21 2019.12.11_raw.txt',
            'Converted SMR 21 2019.12.11.txt',
            n.chamber=3,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR 23 2019.12.09_raw.txt',
            'Converted SMR 23 2019.12.09.txt',
            n.chamber=3,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR 25 2019.12.10_raw.txt',
            'Converted SMR 25 2019.12.10.txt',
            n.chamber=3,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR 25 2019.12.11 take 2_raw.txt',
            'Converted SMR 25 2019.12.11 take 2.txt',
            n.chamber=3,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
#Controls First Run of Upper Penninsula
convert.rMR('empty 2019.12.07_raw.txt', 
            'Converted empty 2019.12.07.txt',
            n.chamber=3,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty 2019.12.12_raw.txt',
            'Converted empty 2019.12.12.txt',
            n.chamber=3,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")


#Import Background Respiration Files
UP1.pre<-import.test('Converted empty 2019.12.07.txt', 
                 info.data=info_08, 
                 logger= "AutoResp",
                 n.chamber=4,
                 plot.temperature= F,
                 plot.oxygen = F)
UP1.post<-import.test('Converted empty 2019.12.12.txt',
                  info.data=info_08, 
                  logger="AutoResp",
                  n.chamber=4,
                  plot.temperature = F,
                  plot.oxygen = F)
#Import Raw SMR Data
UP1.SMR.21.raw<-import.meas('Converted SMR 21 2019.12.08.txt', 
                        info.data=info_08, 
                        logger="AutoResp", 
                        n.chamber=3,
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)
UP1.SMR.21.2.raw<-import.meas('Converted SMR 21 2019.12.11.txt', 
                        info.data=info_08, 
                        n.chamber=3,
                        logger="AutoResp", 
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)
UP1.SMR.23.raw<-import.meas('Converted SMR 23 2019.12.09.txt', 
                          info.data=info_08, 
                          logger="AutoResp", 
                          n.chamber=3,
                          date.format="MDY",
                          plot.temperature = F,
                          plot.oxygen = F)
UP1.SMR.25.raw<-import.meas('Converted SMR 25 2019.12.10.txt', 
                        info.data=info_08, 
                        logger="AutoResp", 
                        n.chamber=3,
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)
UP1.SMR.25.2.raw<-import.meas('Converted SMR 25 2019.12.11 take 2.txt', 
                             info.data=info_08, 
                             logger="AutoResp", 
                             n.chamber=3,
                             date.format="MDY",
                             plot.temperature = F,
                             plot.oxygen = F)
#Clean Raw SMR data with background respiration measurments
UP1.SMR.21.clean<-correct.meas(info.data=info_08,
                           pre.data=UP1.pre,
                           post.data=UP1.post, 
                           meas.data=UP1.SMR.21.raw,
                           method="exponential")
UP1.SMR.21.2.clean<-correct.meas(info.data=info_08,
                           pre.data=UP1.pre,
                           post.data=UP1.post, 
                           meas.data=UP1.SMR.21.2.raw,
                           method="exponential")
UP1.SMR.23.clean<-correct.meas(info.data=info_08,
                             pre.data=UP1.pre,
                             post.data=UP1.post, 
                             meas.data=UP1.SMR.23.raw,
                             method="exponential")
UP1.SMR.25.clean<-correct.meas(info.data=info_08,
                           pre.data=UP1.pre,
                           post.data=UP1.post, 
                           meas.data=UP1.SMR.25.raw,
                           method="exponential")
UP1.SMR.25.2.clean<-correct.meas(info.data=info_08,
                               pre.data=UP1.pre,
                               post.data=UP1.post, 
                               meas.data=UP1.SMR.25.2.raw,
                               method="exponential")
#Extract slopes
UP1.SMR.21.slopes<-extract.slope(UP1.SMR.21.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.25)
UP1.SMR.21.2.slopes<-extract.slope(UP1.SMR.21.2.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.25)
UP1.SMR.23.slopes<-extract.slope(UP1.SMR.23.clean,
                               method="calcSMR.quant",
                               r2=0.95,
                               p=0.25)
UP1.SMR.25.slopes<-extract.slope(UP1.SMR.25.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.25)
UP1.SMR.25.2.slopes<-extract.slope(UP1.SMR.25.2.clean,
                                 method="calcSMR.quant",
                                 r2=0.95,
                                 p=0.25)
#Calculate Metabolic Rates
UP1.SMR.21<-calculate.MR(UP1.SMR.21.slopes,density=1000)
UP1.SMR.21.2<-calculate.MR(UP1.SMR.21.2.slopes,density=1000)
UP1.SMR.23<-calculate.MR(UP1.SMR.23.slopes,density=1000)
UP1.SMR.25<-calculate.MR(UP1.SMR.25.slopes,density=1000)
UP1.SMR.25.2<-calculate.MR(UP1.SMR.25.2.slopes,density=1000)

#Combine Data into one dataframe
UP1.Accute<-rbind.data.frame(UP1.SMR.21,UP1.SMR.21.2,UP1.SMR.23,UP1.SMR.25)
Temperature<-c(21,21,21,21,21,21,23,23,23,25,25,25)
UP1.Accute<-cbind(UP1.Accute, Temperature)
setwd("~/Desktop/Respirometry R code")
write.csv(UP1.Accute,"UP Accute SMR EXP 1")

plot((UP1.Accute$MR.mass/1000)~UP1.Accute$Temperature,
     yaxt="none",
     ylab="Mass specific Metabolic Rate g/L/h",
     xlab= "Temperaturec(C)",
     main="Standard Metabolic Rates", 
     ylim=c(0,.5),
     xlim=c(21,25),
     cex=.75, col=1, pch=19)
axis(2,seq(0,2,.1))





