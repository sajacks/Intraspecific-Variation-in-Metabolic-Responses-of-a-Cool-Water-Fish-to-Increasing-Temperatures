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
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_9.13.2021")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used
info_25.14<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(2.79,2.97,3.27,4.44,2.72,3.07,5.23,5.26), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_25.15<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(2.75,2.83,3.12,4.26,2.59,2.92,4.95,4.93), 
                       Volume = c(78,78,78,78,78,78,78,78), 
                       DO.unit = "mg/L")
info_23.16<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(2.67,2.81,3.04,4.19,2.49,2.83,4.89,4.83), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_21.17<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                                Mass=c(2.64,2.79,2.99,4.12,2.48,2.80,4.83,4.73), 
                                Volume = c(78,78,78,78,78,78,78,78), 
                                DO.unit = "mg/L")
info_25.2.18<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                        Mass=c(2.83,2.84,3.08,4.08,2.47,2.76,4.81,4.65), 
                        Volume = c(78,78,78,78,78,78,78,78), 
                        DO.unit = "mg/L")

#Convert blank runs used for background respiration
convert.rMR('empty_pre2_9.13.2021_raw.txt', 
            'Converted Empty_Pre 09.13.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_9.18.2021_raw.txt',
            'Converted Empty_Post 09.18.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")


#Import Background Respiration Files
DB2.pre<-import.test('Converted Empty_Pre 09.13.2021.txt', 
                     info.data=info_25.14, 
                     logger= "AutoResp",
                     n.chamber=8,
                     plot.temperature= F,
                     plot.oxygen = F)
DB2.post<-import.test('Converted Empty_Post 09.18.2021.txt',
                      info.data=info_25.2.18, 
                      logger="AutoResp",
                      n.chamber=8,
                      plot.temperature = F,
                      plot.oxygen = F)

#######################################################################
#Standard Metabolic Rate
#######################################################################
#Convert %Air Saturation to mg/L for SMR measurements and Empty chamber runs
#First Run of Accute SMR Upper Pennisula
convert.rMR('SMR_25_DB_9.14.2021_raw.txt',
            'Converted SMR 25 09.14.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_23_DB_9.15.2021_raw.txt',
            'Converted SMR 23 09.15.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_21_DB_9.16.2021_raw.txt',
            'Converted SMR 21 09.16.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_25.2_DB_9.17.2021_raw.txt',
            'Converted SMR 25.2 09.17.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
#Import Raw SMR Data
DB2.SMR.25.raw<-import.meas('Converted SMR 25 09.14.2021.txt', 
                        info.data=info_25.14, 
                        logger="AutoResp", 
                        n.chamber=8,
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)
DB2.SMR.23.raw<-import.meas('Converted SMR 23 09.15.2021.txt', 
                          info.data=info_25.15, 
                          logger="AutoResp", 
                          n.chamber=8,
                          date.format="MDY",
                          plot.temperature = F,
                          plot.oxygen = F)
DB2.SMR.21.raw<-import.meas('Converted SMR 21 09.16.2021.txt',
                            info.data=info_23.16,
                            logger="AutoResp",
                            n.chamber=8,date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
DB2.SMR.25.2.raw<-import.meas('Converted SMR 25.2 09.17.2021.txt', 
                              info.data=info_21.17, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)

#Clean Raw SMR data with background respiration measurments
DB2.SMR.25.clean<-correct.meas(info.data=info_25.14,
                               pre.data=DB2.pre,
                               post.data=DB2.post,
                           meas.data=DB2.SMR.25.raw,
                           method="exponential")
DB2.SMR.23.clean<-correct.meas(info.data=info_25.15,
                               pre.data=DB2.pre,
                               post.data=DB2.post, 
                             meas.data=DB2.SMR.23.raw,
                             method="exponential")
DB2.SMR.21.clean<-correct.meas(info.data=info_23.16,
                               pre.data=DB2.pre,
                               post.data=DB2.post,
                           meas.data=DB2.SMR.21.raw,
                           method="exponential")
DB2.SMR.25.2.clean<-correct.meas(info.data=info_21.17,
                                 pre.data=DB2.pre,
                                 post.data=DB2.post, 
                                 meas.data=DB2.SMR.25.2.raw,
                                 method="exponential")

#Extract slopes for Standard Metabolic Rate
DB2.SMR.25.slopes<-extract.slope(DB2.SMR.25.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)

DB2.SMR.23.slopes<-extract.slope(DB2.SMR.23.clean,
                               method="calcSMR.quant",
                               r2=0.95,
                               p=0.2)
DB2.SMR.21.slopes<-extract.slope(DB2.SMR.21.clean,
                             method="calcSMR.quant",
                             r2=0.9,
                             p=0.2)
DB2.SMR.25.2.slopes<-extract.slope(DB2.SMR.25.2.clean,
                                   method="calcSMR.quant",
                                   r2=0.95,
                                   p=0.2)

##Extract slopes for Maximum Metabolic Rate during SMR Measurements
DB2.SMR.MMR.25.slopes<-extract.slope(DB2.SMR.25.clean,
                                 method="max",
                                 n.slope=1)
DB2.SMR.MMR.23.slopes<-extract.slope(DB2.SMR.23.clean,
                                 method="max",
                                 n.slope=1)
DB2.SMR.MMR.21.slopes<-extract.slope(DB2.SMR.21.clean,
                                 method="max",
                                 n.slope=1)
DB2.SMR.MMR.25.2.slopes<-extract.slope(DB2.SMR.25.2.clean,
                                       method="max",
                                       n.slope=1)



##Calculate Standard Metabolic Rates
DB2.SMR.25<-calculate.MR(DB2.SMR.25.slopes,density=1000)
DB2.SMR.23<-calculate.MR(DB2.SMR.23.slopes,density=1000)
DB2.SMR.21<-calculate.MR(DB2.SMR.21.slopes,density=1000)
DB2.SMR.25.2<-calculate.MR(DB2.SMR.25.2.slopes,density=1000)

##Calculate Maximum Metabolic Rates 
DB2.SMR.MMR.25<-calculate.MR(DB2.SMR.MMR.25.slopes,density=1000)
DB2.SMR.MMR.23<-calculate.MR(DB2.SMR.MMR.23.slopes,density=1000)
DB2.SMR.MMR.21<-calculate.MR(DB2.SMR.MMR.21.slopes,density=1000)
DB2.SMR.MMR.25.2<-calculate.MR(DB2.SMR.MMR.25.2.slopes,density=1000)


##Combine Data into one Dataframe
DB2.SMR.Acute<-rbind.data.frame(DB2.SMR.25,DB2.SMR.23,DB2.SMR.21,DB2.SMR.25.2)
DB2.SMR.MMR.Acute<-rbind.data.frame(DB2.SMR.MMR.25,DB2.SMR.MMR.23,DB2.SMR.MMR.21,DB2.SMR.MMR.25.2)

##Add additional Columns
Temperature<-c(25,25,25,25,25,25,25,25,23,23,23,23,23,23,23,23,21,21,21,21,21,21,21,21,25,25,25,25,25,25,25,25)
Trial<-c("SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR","SMR","SMR")
Event<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4)

DB2.SMR.Acute<-cbind(DB2.SMR.Acute,Trial,Temperature,Event)
DB2.SMR.MMR.Acute<-cbind(DB2.SMR.MMR.Acute,Trial,Temperature,Event)

##############33
##Save out CSVs

write.csv(DB2.SMR.Acute,"LP DB Acute SMR EXP 2.csv")
write.csv(DB2.SMR.MMR.Acute,"LP DB Acute SMR-MMR EXP 2.csv")


#######################################################################
#Maximum Metabolic Rate
########################################################################

#Convert %Air Saturation to mg/L for MMR measurements and Empty chamber runs
#First Run of Acute MMR Upper Pennisula
convert.rMR('MMR_25_DB_9.15.2021_raw.txt',
            'Converted MMR 25 09.15.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_23_DB_9.16.2021_raw.txt',
            'Converted MMR 23 09.16.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_21_DB_9.17.2021_raw.txt',
            'Converted MMR 21 09.17.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_25.2_DB_9.18.2021_raw.txt',
            'Converted MMR 25.2 09.18.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw MMR Data
DB2.MMR.25.raw<-import.meas('Converted MMR 25 09.15.2021.txt', 
                            info.data=info_25.15, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
DB2.MMR.23.raw<-import.meas('Converted MMR 23 09.16.2021.txt', 
                            info.data=info_23.16, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
DB2.MMR.21.raw<-import.meas('Converted MMR 21 09.17.2021.txt', 
                            info.data=info_21.17, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
DB2.MMR.25.2.raw<-import.meas('Converted MMR 25.2 09.18.2021.txt', 
                              info.data=info_25.2.18, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)

#Clean Raw MMR data with background respiration measurments
DB2.MMR.25.clean<-correct.meas(info.data=info_25.15,
                               pre.data=DB2.pre,
                               post.data=DB2.post,
                               meas.data=DB2.MMR.25.raw,
                               method="exponential")
                                
DB2.MMR.23.clean<-correct.meas(info.data=info_23.16,
                               pre.data=DB2.pre,
                               post.data=DB2.post, 
                               meas.data=DB2.MMR.23.raw,
                               method="exponential")
DB2.MMR.21.clean<-correct.meas(info.data=info_21.17,
                               pre.data=DB2.pre,
                               post.data=DB2.post,
                               meas.data=DB2.MMR.21.raw,
                               method="exponential")
DB2.MMR.25.2.clean<-correct.meas(info.data=info_25.2.18,
                                 pre.data=DB2.pre,
                                 post.data=DB2.post, 
                                 meas.data=DB2.MMR.25.2.raw,
                                 method="exponential")

#Extract slopes
DB2.MMR.25.slopes<-extract.slope(DB2.MMR.25.clean,
                                 method="max",
                                 n.slope=1)
DB2.MMR.23.slopes<-extract.slope(DB2.MMR.23.clean,
                                 method="max",
                                 n.slope=1)
DB2.MMR.21.slopes<-extract.slope(DB2.MMR.21.clean,
                                 method="max",
                                 n.slope=1)
DB2.MMR.25.2.slopes<-extract.slope(DB2.MMR.25.2.clean,
                                   method="max",
                                   n.slope=1)

##Calculate Metabolic Rates
DB2.MMR.25<-calculate.MR(DB2.MMR.25.slopes,density=1000)
DB2.MMR.23<-calculate.MR(DB2.MMR.23.slopes,density=1000)
DB2.MMR.21<-calculate.MR(DB2.MMR.21.slopes,density=1000)
DB2.MMR.25.2<-calculate.MR(DB2.MMR.25.2.slopes,density=1000)


##Combine data into one Dataframe
DB2.MMR.Acute<-rbind.data.frame(DB2.MMR.25,DB2.MMR.23,DB2.MMR.21,DB2.MMR.25.2)

##Add additional columns
Trial<-c("MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR","MMR","MMR")
DB2.MMR.Acute<-cbind(DB2.MMR.Acute,Trial,Temperature,Event)

##Save data our as CSV

write.csv(DB2.MMR.Acute,"LP DB Acute MMR EXP 2.csv")

#######################################################################
#Data Analysis and grapahs
#######################################################################

DB2.MMR<-rbind(DB2.SMR.MMR.Acute,DB2.MMR.Acute)

DB2.MMR<-read.csv("LP DB Acute MMR BEST EXP 2.csv")
DB2.MMR$Temperature<-factor(Temperature)
DB2.SMR.MMR.Acute$Temperature<-factor(Temperature)
DB2.MMR.Acute$Temperature<-factor(Temperature)
MMRC$Temperature<-factor(Temperature)
MMRC<-as.data.frame(MMRC)
m<-ggplot() 
m<-m+geom_point(data=DB2.SMR.MMR.Acute,aes(Temperature, MR.mass/60, group=Ind,color=Ind)) + ggtitle("DB4 MMR, Comparison")+theme_classic()
m<-m+geom_point(data=DB2.MMR.Acute,aes(Temperature, MR.mass/60, group=Ind,color=Ind))
m<-m+ylab("MO2 (mgO2/(kg min))")
m

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
        ggplot(DB2.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
                geom_boxplot() +geom_boxplot(data=DB2.MMR, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
                ylab("MO2 (mgO2/(kg min))") + ggtitle("DB SMR,Best MMR ") + theme(legend.position="none"),
ggplot(DB2.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
        geom_boxplot() +geom_boxplot(data=DB2.MMR.Acute, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
        ylab("MO2 (mgO2/(kg min))") + ggtitle("DB SMR, MMR ") + theme(legend.position="none"),
ggplot(DB2.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
        geom_boxplot() +geom_boxplot(data=DB2.SMR.MMR.Acute, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
        ylab("MO2 (mgO2/(kg min))") + ggtitle("DB SMR, SMR-MMR ") + theme(legend.position="none"),
        
        ncol=3
)

grid.arrange(
        ggplot(subset(DB2.MMR,Event%in%c(1)), 
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


###################################################################
####Additional Analysis for Manuscript
###################################################################

###Function from Chabot et al 2016 appendix 1

calcSMR = function(Y, q=c(0.1,0.15,0.2,0.25,0.3), G=1:4){
  u = sort(Y)
  the.Mclust <- Mclust(Y,  G=G)
  cl <- the.Mclust$classification
  # sometimes, the class containing SMR is not called 1
  # the following presumes that when class 1 contains > 10% of cases, 
  # it contains SMR, otherwise we take class 2
  cl2 <- as.data.frame(table(cl))
  cl2$cl <- as.numeric(levels(cl2$cl))
  valid <- cl2$Freq>=0.1*length(time)  
  the.cl <- min(cl2$cl[valid])
  left.distr <- Y[the.Mclust$classification==the.cl]
  mlnd = the.Mclust$parameters$mean[the.cl]
  CVmlnd = sd(left.distr)/mlnd * 100
  quant=quantile(Y, q)
  low10=mean(u[1:10])
  low10pc = mean(u[6:(5 + round(0.1*(length(u)-5)))])
  # remove 5 outliers, keep lowest 10% of the rest, average
  # Herrmann & Enders 2000
  return(list(mlnd=mlnd, quant=quant, low10=low10, low10pc=low10pc,
              cl=cl, CVmlnd=CVmlnd))
}






##Extract all slopes for Standard Metabolic Rate
DB2.SMR.21.all<-extract.slope(DB2.SMR.21.clean,
                              method="all",
                              r2=0.9)
##Calculate MR for all slopes
DB2.SMR.21<-calculate.MR(DB2.SMR.21.all,density=1000)
##Turn Phase into a numeric variable
DB2.SMR.21<-DB2.SMR.21%>%mutate(Phase=as.numeric(gsub("M","",Phase)))

##Extract average phase used for calculating SMR
grouped_quantiles<-DB2.SMR.21 %>%
  group_by(Ind) %>%
  summarise(twentieth_quantile =quantile(MR.abs,.2,na.rm=TRUE)) %>%
  ungroup()

extract_values<-function(group_df,q){
  group_df %>%
    filter(MR.abs<= q)
}

values<-list()

for (i in seq_along(grouped_quantiles$Ind)){
  group_data_subset<-DB2.SMR.21%>%
    filter(Ind==grouped_quantiles$Ind[i])
  quantile_value<-grouped_quantiles$twentieth_quantile[i]
  values[[grouped_quantiles$Ind[i]]]<-extract_values(group_data_subset,quantile_value)
}

for( i in 1:length(values)) {
  group_df<-values[[i]]
  mean_phase_values[i]<-mean(group_df$Phase,na.rm=TRUE)
}
DB2.21<-mean_phase_values



DB2.SMR.23.all<-extract.slope(DB2.SMR.23.clean,
                              method="all",
                              r2=0.9)
DB2.SMR.23<-calculate.MR(DB2.SMR.23.all,density=1000)
DB2.SMR.23<-DB2.SMR.23%>%mutate(Phase=as.numeric(gsub("M","",Phase)))

##Extract average phase used for calculating SMR
grouped_quantiles<-DB2.SMR.23 %>%
  group_by(Ind) %>%
  summarise(twentieth_quantile =quantile(MR.abs,.2,na.rm=TRUE)) %>%
  ungroup()

extract_values<-function(group_df,q){
  group_df %>%
    filter(MR.abs<= q)
}

values<-list()

for (i in seq_along(grouped_quantiles$Ind)){
  group_data_subset<-DB2.SMR.23%>%
    filter(Ind==grouped_quantiles$Ind[i])
  quantile_value<-grouped_quantiles$twentieth_quantile[i]
  values[[grouped_quantiles$Ind[i]]]<-extract_values(group_data_subset,quantile_value)
}

for( i in 1:length(values)) {
  group_df<-values[[i]]
  mean_phase_values[i]<-mean(group_df$Phase,na.rm=TRUE)
}
DB2.23<-mean_phase_values


DB2.SMR.25.all<-extract.slope(DB2.SMR.25.clean,
                              method="all",
                              r2=0.9)
DB2.SMR.25<-calculate.MR(DB2.SMR.25.all,density=1000)
DB2.SMR.25<-DB2.SMR.25%>%mutate(Phase=as.numeric(gsub("M","",Phase)))

##Extract average phase used for calculating SMR
grouped_quantiles<-DB2.SMR.25 %>%
  group_by(Ind) %>%
  summarise(twentieth_quantile =quantile(MR.abs,.2,na.rm=TRUE)) %>%
  ungroup()

extract_values<-function(group_df,q){
  group_df %>%
    filter(MR.abs<= q)
}

values<-list()

for (i in seq_along(grouped_quantiles$Ind)){
  group_data_subset<-DB2.SMR.25%>%
    filter(Ind==grouped_quantiles$Ind[i])
  quantile_value<-grouped_quantiles$twentieth_quantile[i]
  values[[grouped_quantiles$Ind[i]]]<-extract_values(group_data_subset,quantile_value)
}

for( i in 1:length(values)) {
  group_df<-values[[i]]
  mean_phase_values[i]<-mean(group_df$Phase,na.rm=TRUE)
}
DB2.25<-mean_phase_values


Phase_values<-data.frame(rbind(DB2.21,DB2.23,DB2.25),row.names=c(1,2,3))

Phase_values<-cbind(Phase_values,Temperature)
names(Phase_values)<-c(1,2,3,4,5,6,7,8,"Temp")
library(reshape2)
Phase<-melt(Phase_values,id.vars="Temp")
ggplot(Phase,aes(x=variable,y=value,color=Temp))+geom_point()

install.packages("quantreg")
library("quantreg")
summary(rq(MR.abs~Phase, tau=.2,data=DB2.SMR.21),se="nid")

summary(rq(MR.abs~poly(Phase,2)*Ind, tau=.2,data=DB2.SMR.23),se="nid")

summary(rq(MR.abs~poly(Phase,2)*Ind, tau=.2,data=DB2.SMR.25),se="nid")

summary(rq(MR.abs~poly(Phase,2)*Ind, tau=.2,data=DB2.SMR.21.2),se="nid")

install.packages("lqmm")
DB2.SMR.21$Ind<-as.numeric(DB2.SMR.21$Ind)
DB2.SMR.23$Ind<-as.numeric(DB2.SMR.23$Ind)
DB2.SMR.25$Ind<-as.numeric(DB2.SMR.25$Ind)
DB2.SMR.25.2$Ind<-as.numeric(DB2.SMR.25.2$Ind)
library(lqmm)
DB2.21.model<-lqmm(fixed=MR.abs~Phase,random=~1|Ind,group=Ind,tau=.2,data=DB2.SMR.21)
summary(DB2.21.model)

DB2.23.model<-lqmm(fixed=MR.abs~Phase,random=~1|Ind,group=Ind,tau=.2,data=DB2.SMR.23)
summary(DB2.23.model)

DB2.25.model<-lqmm(fixed=MR.abs~Phase,random=~1|Ind,group=Ind,tau=.2,data=DB2.SMR.25,control = list(eps = 1e-2,  # Increase tolerance by setting a higher eps value
                                                                                                    LP_max_iter = 1000))
summary(DB2.25.model)

DB2.25.2.model<-lqmm(fixed=MR.abs~Phase,random=~1|Ind,group=Ind,tau=.2,data=DB2.SMR.25.2,control = list(eps = 1e-2,  # Increase tolerance by setting a higher eps value
                                                                                                    LP_max_iter = 1000,verbose=T))
summary(DB2.25.2.model)

DB2.SMR.25.2.all<-extract.slope(DB2.SMR.25.2.clean,
                                method="all",
                                r2=0.9)
DB2.SMR.25.2<-calculate.MR(DB2.SMR.25.2.all,density=1000)
DB2.SMR.25.2<-DB2.SMR.25.2%>%mutate(Phase=as.numeric(gsub("M","",Phase)))


ggplot(DB2.SMR.21,aes(x=Phase, y=MR.abs,color=Chamber.No))+geom_point()+geom_smooth(method=lm,formula=y~poly(x,2))
grid.arrange(
  ggplot(DB2.SMR.21,aes(x=Phase, y=MR.abs,color=Chamber.No))+geom_point()+geom_smooth(method=lm,formula=y~poly(x,2)),
  ggplot(DB2.SMR.23,aes(x=Phase, y=MR.abs,color=Chamber.No))+geom_point()+geom_smooth(method=lm,formula=y~poly(x,2)),
  ggplot(DB2.SMR.25,aes(x=Phase, y=MR.abs,color=Chamber.No))+geom_point()+geom_smooth(method=lm,formula=y~poly(x,2)),
  ggplot(DB2.SMR.21.2,aes(x=Phase, y=MR.abs,color=Chamber.No))+geom_point()+geom_smooth(method=lm,formula=y~poly(x,2)),
  ncol=2
)

write.csv(DB2.SMR.21,"DB2.SMR.21.all.csv")
write.csv(DB2.SMR.23,"DB2.SMR.23.all.csv")
write.csv(DB2.SMR.25,"DB2.SMR.25.all.csv")
write.csv(DB2.SMR.21.2,"DB2.SMR.21.2.all.csv")




