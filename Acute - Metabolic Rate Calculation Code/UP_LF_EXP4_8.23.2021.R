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
setwd("~/Desktop/MR data/2021 Spring Fingerling Data/Static MR_8.23.2021")
#Input Chamber ID information, Fish Wet weight, Volume of chambers and Tubing, Final DO unit used
info_25.23<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(5.15,4.66,5.11,4.92,5.40,5.41,3.27,5.14), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_25.24<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                       Mass=c(4.99,4.35,4.97,4.82,5.18,5.37,3.12,4.98), 
                       Volume = c(78,78,78,78,78,78,78,78), 
                       DO.unit = "mg/L")
info_23.25<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                    Mass=c(4.94,4.4,4.95,4.66,5.01,5.14,3.01,4.81), 
                    Volume = c(78,78,78,78,78,78,78,78), 
                    DO.unit = "mg/L")
info_21.26<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                                Mass=c(5.06,4.25,4.99,4.58,4.94,4.99,3.01,4.76), 
                                Volume = c(78,78,78,78,78,78,78,78), 
                                DO.unit = "mg/L")
info_25.2.27<-input.info(ID=c(1,2,3,4,5,6,7,8), 
                        Mass=c(4.83,4.18,4.72,4.66,4.8,4.95,2.98,4.73), 
                        Volume = c(78,78,78,78,78,78,78,78), 
                        DO.unit = "mg/L")

#Convert blank runs used for background respiration
convert.rMR('empty_pre_8.23.2021_raw.txt', 
            'Converted Empty_Pre 08.23.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('empty_post_8.27.2021_raw.txt',
            'Converted Empty_Post 08.27.2021.txt',
            n.chamber=8,
            logger="AutoResp",
            DO.units.in="pct",
            DO.units.out="mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")


#Import Background Respiration Files
LF4.pre<-import.test('Converted Empty_Pre 08.23.2021.txt', 
                     info.data=info_21.24, 
                     logger= "AutoResp",
                     n.chamber=8,
                     plot.temperature= F,
                     plot.oxygen = F)
LF4.post<-import.test('Converted Empty_Post 08.27.2021.txt',
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
convert.rMR('SMR_25_LF_8.23.2021_raw.txt',
            'Converted SMR 25 08.23.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_23_LF_8.24.2021_raw.txt',
            'Converted SMR 23 08.24.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_21_LF_8.25.2021_raw.txt',
            'Converted SMR 21 08.25.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('SMR_25.2_LF_8.26.2021_raw.txt',
            'Converted SMR 25.2 08.26.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
#Import Raw SMR Data
LF4.SMR.25.raw<-import.meas('Converted SMR 25 08.23.2021.txt', 
                        info.data=info_21.24, 
                        logger="AutoResp", 
                        n.chamber=8,
                        date.format="MDY",
                        plot.temperature = F,
                        plot.oxygen = F)
LF4.SMR.23.raw<-import.meas('Converted SMR 23 08.24.2021.txt', 
                          info.data=info_23.25, 
                          logger="AutoResp", 
                          n.chamber=8,
                          date.format="MDY",
                          plot.temperature = F,
                          plot.oxygen = F)
LF4.SMR.21.raw<-import.meas('Converted SMR 21 08.25.2021.txt',
                            info.data=info_25.26,
                            logger="AutoResp",
                            n.chamber=8,date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
LF4.SMR.25.2.raw<-import.meas('Converted SMR 25.2 08.26.2021.txt', 
                              info.data=info_21.2.27, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)

#Clean Raw SMR data with background respiration measurments
LF4.SMR.25.clean<-correct.meas(info.data=info_25.24,
                               pre.data=LF4.pre,
                               post.data=LF4.post,
                           meas.data=LF4.SMR.25.raw,
                           method="exponential")
LF4.SMR.23.clean<-correct.meas(info.data=info_23.25,
                               pre.data=LF4.pre,
                               post.data=LF4.post, 
                             meas.data=LF4.SMR.23.raw,
                             method="exponential")
LF4.SMR.21.clean<-correct.meas(info.data=info_21.26,
                               pre.data=LF4.pre,
                               post.data=LF4.post,
                           meas.data=LF4.SMR.21.raw,
                           method="exponential")
LF4.SMR.25.2.clean<-correct.meas(info.data=info_25.2.27,
                                 pre.data=LF4.pre,
                                 post.data=LF4.post, 
                                 meas.data=LF4.SMR.25.2.raw,
                                 method="exponential")

#Extract slopes for Standard Metabolic Rate
LF4.SMR.25.slopes<-extract.slope(LF4.SMR.25.clean,
                             method="calcSMR.quant",
                             r2=0.95,
                             p=0.2)

LF4.SMR.23.slopes<-extract.slope(LF4.SMR.23.clean,
                               method="calcSMR.quant",
                               r2=0.95,
                               p=0.2)
LF4.SMR.21.slopes<-extract.slope(LF4.SMR.21.clean,
                             method="calcSMR.quant",
                             r2=0.9,
                             p=0.2)
LF4.SMR.25.2.slopes<-extract.slope(LF4.SMR.25.2.clean,
                                   method="calcSMR.quant",
                                   r2=0.95,
                                   p=0.2)

##Extract slopes for Maximum Metabolic Rate during SMR Measurements
LF4.SMR.MMR.25.slopes<-extract.slope(LF4.SMR.25.clean,
                                 method="max",
                                 n.slope=1)
LF4.SMR.MMR.23.slopes<-extract.slope(LF4.SMR.23.clean,
                                 method="max",
                                 n.slope=1)
LF4.SMR.MMR.21.slopes<-extract.slope(LF4.SMR.21.clean,
                                 method="max",
                                 n.slope=1)
LF4.SMR.MMR.25.2.slopes<-extract.slope(LF4.SMR.25.2.clean,
                                       method="max",
                                       n.slope=1)



##Calculate Standard Metabolic Rates
LF4.SMR.25<-calculate.MR(LF4.SMR.25.slopes,density=1000)
LF4.SMR.23<-calculate.MR(LF4.SMR.23.slopes,density=1000)
LF4.SMR.21<-calculate.MR(LF4.SMR.21.slopes,density=1000)
LF4.SMR.25.2<-calculate.MR(LF4.SMR.25.2.slopes,density=1000)

##Calculate Maximum Metabolic Rates 
LF4.SMR.MMR.25<-calculate.MR(LF4.SMR.MMR.25.slopes,density=1000)
LF4.SMR.MMR.23<-calculate.MR(LF4.SMR.MMR.23.slopes,density=1000)
LF4.SMR.MMR.21<-calculate.MR(LF4.SMR.MMR.21.slopes,density=1000)
LF4.SMR.MMR.25.2<-calculate.MR(LF4.SMR.MMR.25.2.slopes,density=1000)


##Combine Data into one Dataframe
LF4.SMR.Acute<-rbind.data.frame(LF4.SMR.25,LF4.SMR.23,LF4.SMR.21,LF4.SMR.25.2)
LF4.SMR.MMR.Acute<-rbind.data.frame(LF4.SMR.MMR.25,LF4.SMR.MMR.23,LF4.SMR.MMR.21,LF4.SMR.MMR.25.2)

##Add additional Columns
Temperature<-c(25,25,25,25,25,25,25,25,23,23,23,23,23,23,23,23,21,21,21,21,21,21,21,21,25,25,25,25,25,25,25,25)
Trial<-c("SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR",
         "SMR","SMR","SMR","SMR","SMR","SMR","SMR")
Event<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4)

LF4.SMR.Acute<-cbind(LF4.SMR.Acute,Trial,Temperature,Event)
LF4.SMR.MMR.Acute<-cbind(LF4.SMR.MMR.Acute,Trial,Temperature,Event)

##############33
##Save out CSVs

write.csv(LF4.SMR.Acute,"UP LF Acute SMR EXP 4.csv")
write.csv(LF4.SMR.MMR.Acute,"UP LF Acute SMR-MMR EXP 4.csv")


#######################################################################
#Maximum Metabolic Rate
########################################################################

#Convert %Air Saturation to mg/L for MMR measurements and Empty chamber runs
#First Run of Acute MMR Upper Pennisula
convert.rMR('MMR_25_LF_8.24.2021_raw.txt',
            'Converted MMR 25 08.24.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_23_LF_8.25.2021_raw.txt',
            'Converted MMR 23 08.25.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_21_LF_8.26.2021_raw.txt',
            'Converted MMR 21 08.26.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")
convert.rMR('MMR_25.2_LF_8.27.2021_raw.txt',
            'Converted MMR 25.2 08.27.2021.txt',
            n.chamber=8,logger="AutoResp",
            DO.units.in = "pct",
            DO.units.out = "mg/L",
            salinity=0,
            bar.press=101.3,
            bar.units.in="kpa")

#Import Raw MMR Data
LF4.MMR.25.raw<-import.meas('Converted MMR 25 08.24.2021.txt', 
                            info.data=info_21.24, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
LF4.MMR.23.raw<-import.meas('Converted MMR 23 08.25.2021.txt', 
                            info.data=info_23.25, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
LF4.MMR.21.raw<-import.meas('Converted MMR 21 08.26.2021.txt', 
                            info.data=info_21.26, 
                            logger="AutoResp", 
                            n.chamber=8,
                            date.format="MDY",
                            plot.temperature = F,
                            plot.oxygen = F)
LF4.MMR.25.2.raw<-import.meas('Converted MMR 25.2 08.27.2021.txt', 
                              info.data=info_25.2.27, 
                              n.chamber=8,
                              logger="AutoResp", 
                              date.format="MDY",
                              plot.temperature = F,
                              plot.oxygen = F)

#Clean Raw MMR data with background respiration measurments
LF4.MMR.25.clean<-correct.meas(info.data=info_25.24,
                               pre.data=LF4.pre,
                               post.data=LF4.post,
                               meas.data=LF4.MMR.25.raw,
                               method="exponential")
                                
LF4.MMR.23.clean<-correct.meas(info.data=info_23.25,
                               pre.data=LF4.pre,
                               post.data=LF4.post, 
                               meas.data=LF4.MMR.23.raw,
                               method="exponential")
LF4.MMR.21.clean<-correct.meas(info.data=info_21.26,
                               pre.data=LF4.pre,
                               post.data=LF4.post,
                               meas.data=LF4.MMR.21.raw,
                               method="exponential")
LF4.MMR.25.2.clean<-correct.meas(info.data=info_25.2.27,
                                 pre.data=LF4.pre,
                                 post.data=LF4.post, 
                                 meas.data=LF4.MMR.25.2.raw,
                                 method="exponential")

#Extract slopes
LF4.MMR.25.slopes<-extract.slope(LF4.MMR.25.clean,
                                 method="max",
                                 n.slope=1)
LF4.MMR.23.slopes<-extract.slope(LF4.MMR.23.clean,
                                 method="max",
                                 n.slope=1)
LF4.MMR.21.slopes<-extract.slope(LF4.MMR.21.clean,
                                 method="max",
                                 n.slope=1)
LF4.MMR.25.2.slopes<-extract.slope(LF4.MMR.25.2.clean,
                                   method="max",
                                   n.slope=1)

##Calculate Metabolic Rates
LF4.MMR.25<-calculate.MR(LF4.MMR.25.slopes,density=1000)
LF4.MMR.23<-calculate.MR(LF4.MMR.23.slopes,density=1000)
LF4.MMR.21<-calculate.MR(LF4.MMR.21.slopes,density=1000)
LF4.MMR.25.2<-calculate.MR(LF4.MMR.25.2.slopes,density=1000)


##Combine data into one Dataframe
LF4.MMR.Acute<-rbind.data.frame(LF4.MMR.25,LF4.MMR.23,LF4.MMR.21,LF4.MMR.25.2)

##Add additional columns
Trial<-c("MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR",
         "MMR","MMR","MMR","MMR","MMR","MMR","MMR")
LF4.MMR.Acute<-cbind(LF4.MMR.Acute,Trial,Temperature,Event)

##Save data our as CSV
setwd("~/Desktop/Dissertation/Chapter 2")

write.csv(LF4.MMR.Acute,"UP LF Acute MMR EXP 4.csv")

#######################################################################
#Data Analysis and grapahs
#######################################################################

LF4.MMR<-rbind(LF4.SMR.MMR.Acute,LF4.MMR.Acute)

LF4.MMR<-read.csv("UP LF Acute MMR BEST EXP 3.csv")
LF4.MMR$Temperature<-factor(Temperature)
LF4.SMR.MMR.Acute$Temperature<-factor(Temperature)
LF4.MMR.Acute$Temperature<-factor(Temperature)
MMRC$Temperature<-factor(Temperature)
MMRC<-as.data.frame(MMRC)
m<-ggplot() 
m<-m+geom_point(data=LF4.SMR.MMR.Acute,aes(Temperature, MR.mass/60, group=Ind,color=Ind)) + ggtitle("LF4 MMR, Comparison")+theme_classic()
m<-m+geom_point(data=LF4.MMR.Acute,aes(Temperature, MR.mass/60, group=Ind,color=Ind))
m<-m+ylab("MO2 (mgO2/(kg min))")
m

s<-ggplot(LF2.SMR,aes(Trial, MR.mass/60, group=Ind,color=Ind)) 
s<-s+geom_line() + geom_point() + ggtitle("LF2 SMR, Comparison")+theme_classic()
s<-s+ylab("MO2 (mgO2/(kg min))")
s
LF4.SMR.Acute$Temperature<-factor(Temperature)

SMR<-ggplot(LF4.SMR.Acute,aes(x=Temperature,y=MR.mass/60))
SMR<-SMR+geom_boxplot() + ggtitle("LF4 SMR, SMR")+theme_classic()
SMR<-SMR+geom_boxplot(data=LF4.MMR, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ylab("MO2 (mgO2/(kg min))")
SMR


grid.arrange(
        ggplot(LF4.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
                geom_boxplot() +geom_boxplot(data=LF4.MMR, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
                ylab("MO2 (mgO2/(kg min))") + ggtitle("LF SMR,Best MMR ") + theme(legend.position="none"),
ggplot(LF4.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
        geom_boxplot() +geom_boxplot(data=LF4.MMR.Acute, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
        ylab("MO2 (mgO2/(kg min))") + ggtitle("LF SMR, MMR ") + theme(legend.position="none"),
ggplot(LF4.SMR.Acute,aes(x=Temperature,y=MR.mass/60)) + 
        geom_boxplot() +geom_boxplot(data=LF4.SMR.MMR.Acute, aes(x=as.factor(Temperature),y=MR.mass/60,color="red"))+ 
        ylab("MO2 (mgO2/(kg min))") + ggtitle("LF SMR, SMR-MMR ") + theme(legend.position="none"),
        
        ncol=3
)

grid.arrange(
        ggplot(subset(LF4.MMR,Event%in%c(1)), 
               aes(x = Trial, y = MR.mass, group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF 25") + theme(legend.position="none")+ylim(200,1500),
        ggplot(subset(LF4.MMR,Event%in%c(4)), 
               aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF 25.2") + theme(legend.position="none")+ylim(200,1500),
        ggplot(subset(LF4.MMR,Event%in%c(2)), 
               aes(x = Trial, y = MR.mass, group=Ind, color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF 23") + theme(legend.position="none")+ylim(200,1500),
        ggplot(subset(LF4.MMR,Event%in%c(3)), 
               aes(x = Trial, y = MR.mass,  group=Ind,color=Ind)) + 
                geom_point()+geom_line() + ggtitle("MMR, LF 21") + theme(legend.position="none")+ylim(200,1500),
        
        ncol=2
)
aMMR<-aov(MMRC$MR.mass~MMRC$Temperature*MMRC$Trial,data=MMRC)
Anova(aMMR)
TukeyHSD((aMMR))

as.factor(LF4.MMR$Temperature)


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
LF4.SMR.21.all<-extract.slope(LF4.SMR.21.clean,
                              method="all",
                              r2=0.9)
##Calculate MR for all slopes
LF4.SMR.21<-calculate.MR(LF4.SMR.21.all,density=1000)
##Turn Phase into a numeric variable
LF4.SMR.21<-LF4.SMR.21%>%mutate(Phase=as.numeric(gsub("M","",Phase)))

##Extract average phase used for calculating SMR
grouped_quantiles<-LF4.SMR.21 %>%
  group_by(Ind) %>%
  summarise(twentieth_quantile =quantile(MR.abs,.2,na.rm=TRUE)) %>%
  ungroup()

extract_values<-function(group_df,q){
  group_df %>%
    filter(MR.abs<= q)
}

values<-list()

for (i in seq_along(grouped_quantiles$Ind)){
  group_data_subset<-LF4.SMR.21%>%
    filter(Ind==grouped_quantiles$Ind[i])
  quantile_value<-grouped_quantiles$twentieth_quantile[i]
  values[[grouped_quantiles$Ind[i]]]<-extract_values(group_data_subset,quantile_value)
}

for( i in 1:length(values)) {
  group_df<-values[[i]]
  mean_phase_values[i]<-mean(group_df$Phase,na.rm=TRUE)
}
LF4.21<-mean_phase_values



LF4.SMR.23.all<-extract.slope(LF4.SMR.23.clean,
                              method="all",
                              r2=0.9)
LF4.SMR.23<-calculate.MR(LF4.SMR.23.all,density=1000)
LF4.SMR.23<-LF4.SMR.23%>%mutate(Phase=as.numeric(gsub("M","",Phase)))

##Extract average phase used for calculating SMR
grouped_quantiles<-LF4.SMR.23 %>%
  group_by(Ind) %>%
  summarise(twentieth_quantile =quantile(MR.abs,.2,na.rm=TRUE)) %>%
  ungroup()

extract_values<-function(group_df,q){
  group_df %>%
    filter(MR.abs<= q)
}

values<-list()

for (i in seq_along(grouped_quantiles$Ind)){
  group_data_subset<-LF4.SMR.23%>%
    filter(Ind==grouped_quantiles$Ind[i])
  quantile_value<-grouped_quantiles$twentieth_quantile[i]
  values[[grouped_quantiles$Ind[i]]]<-extract_values(group_data_subset,quantile_value)
}

for( i in 1:length(values)) {
  group_df<-values[[i]]
  mean_phase_values[i]<-mean(group_df$Phase,na.rm=TRUE)
}
LF4.23<-mean_phase_values


LF4.SMR.25.all<-extract.slope(LF4.SMR.25.clean,
                              method="all",
                              r2=0.9)
LF4.SMR.25<-calculate.MR(LF4.SMR.25.all,density=1000)
LF4.SMR.25<-LF4.SMR.25%>%mutate(Phase=as.numeric(gsub("M","",Phase)))

##Extract average phase used for calculating SMR
grouped_quantiles<-LF4.SMR.25 %>%
  group_by(Ind) %>%
  summarise(twentieth_quantile =quantile(MR.abs,.2,na.rm=TRUE)) %>%
  ungroup()

extract_values<-function(group_df,q){
  group_df %>%
    filter(MR.abs<= q)
}

values<-list()

for (i in seq_along(grouped_quantiles$Ind)){
  group_data_subset<-LF4.SMR.25%>%
    filter(Ind==grouped_quantiles$Ind[i])
  quantile_value<-grouped_quantiles$twentieth_quantile[i]
  values[[grouped_quantiles$Ind[i]]]<-extract_values(group_data_subset,quantile_value)
}

for( i in 1:length(values)) {
  group_df<-values[[i]]
  mean_phase_values[i]<-mean(group_df$Phase,na.rm=TRUE)
}
LF4.25<-mean_phase_values


Phase_values<-data.frame(rbind(LF4.21,LF4.23,LF4.25),row.names=c(1,2,3))

Phase_values<-cbind(Phase_values,Temperature)
names(Phase_values)<-c(1,2,3,4,5,6,7,8,"Temp")
library(reshape2)
Phase<-melt(Phase_values,id.vars="Temp")
ggplot(Phase,aes(x=variable,y=value,color=Temp))+geom_point()

install.packages("quantreg")
library("quantreg")
summary(rq(MR.abs~Phase, tau=.2,data=LF4.SMR.21),se="nid")

summary(rq(MR.abs~poly(Phase,2)*Ind, tau=.2,data=LF4.SMR.23),se="nid")

summary(rq(MR.abs~poly(Phase,2)*Ind, tau=.2,data=LF4.SMR.25),se="nid")

summary(rq(MR.abs~poly(Phase,2)*Ind, tau=.2,data=LF4.SMR.21.2),se="nid")

install.packages("lqmm")
LF4.SMR.21$Ind<-as.numeric(LF4.SMR.21$Ind)
LF4.SMR.23$Ind<-as.numeric(LF4.SMR.23$Ind)
LF4.SMR.25$Ind<-as.numeric(LF4.SMR.25$Ind)
LF4.SMR.25.2$Ind<-as.numeric(LF4.SMR.25.2$Ind)
library(lqmm)
LF4.21.model<-lqmm(fixed=MR.abs~Phase,random=~1|Ind,group=Ind,tau=.2,data=LF4.SMR.21)
summary(LF4.21.model)

LF4.23.model<-lqmm(fixed=MR.abs~Phase,random=~1|Ind,group=Ind,tau=.2,data=LF4.SMR.23)
summary(LF4.23.model)

LF4.25.model<-lqmm(fixed=MR.abs~Phase,random=~1|Ind,group=Ind,tau=.2,data=LF4.SMR.25,control = list(eps = 1e-4,  # Increase tolerance by setting a higher eps value
                                                                                                    LP_max_iter = 1000))
summary(LF4.25.model)

LF4.25.2.model<-lqmm(fixed=MR.abs~Phase,random=~1|Ind,group=Ind,tau=.2,data=LF4.SMR.25.2,control = list(eps = 1e-4,  # Increase tolerance by setting a higher eps value
                                                                                                    LP_max_iter = 1000))
summary(LF4.25.2.model)


LF4.SMR.25.2.all<-extract.slope(LF4.SMR.25.2.clean,
                                method="all",
                                r2=0.9)
LF4.SMR.25.2<-calculate.MR(LF4.SMR.25.2.all,density=1000)
LF4.SMR.25.2<-LF4.SMR.25.2%>%mutate(Phase=as.numeric(gsub("M","",Phase)))


ggplot(LF4.SMR.21,aes(x=Phase, y=MR.abs,color=Chamber.No))+geom_point()+geom_smooth(method=lm,formula=y~poly(x,2))
grid.arrange(
  ggplot(LF4.SMR.21,aes(x=Phase, y=MR.abs,color=Chamber.No))+geom_point()+geom_smooth(method=lm,formula=y~poly(x,2)),
  ggplot(LF4.SMR.23,aes(x=Phase, y=MR.abs,color=Chamber.No))+geom_point()+geom_smooth(method=lm,formula=y~poly(x,2)),
  ggplot(LF4.SMR.25,aes(x=Phase, y=MR.abs,color=Chamber.No))+geom_point()+geom_smooth(method=lm,formula=y~poly(x,2)),
  ggplot(LF4.SMR.21.2,aes(x=Phase, y=MR.abs,color=Chamber.No))+geom_point()+geom_smooth(method=lm,formula=y~poly(x,2)),
  ncol=2
)

write.csv(LF4.SMR.21,"LF4.SMR.21.all.csv")
write.csv(LF4.SMR.23,"LF4.SMR.23.all.csv")
write.csv(LF4.SMR.25,"LF4.SMR.25.all.csv")
write.csv(LF4.SMR.21.2,"LF4.SMR.21.2.all.csv")





