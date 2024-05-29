pkgs <- c("here", "readr", "janitor", "mgcv", "gratia", "dplyr", "ggplot2","glmm",
          "ggrepel","visreg","tibble","mgcViz","tidymv","tidyverse", "nlme", "emmeans", "kableExtra", "GGally", 
          "qqplotr","grid","gridExtra","car","RColorBrewer","r2glmm","lmerTest","ggpubr","multcomp","jtools","lme4","sjPlot","MuMIn","itsadug")

vapply(pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)

#Set working directory
setwd("~/Desktop/Chapter 2 Manuscript")

#Read in Raw Data
all_data <- read.csv("Metabolism.Data.csv")
#Check data structure
head(all_data)

sapply(all_data,class)
sapply(final_data,class)

all_data<- transform(
  all_data,
  Ind=as.factor(Ind),
  SMR.units=as.numeric(SMR.units),
  test.Num=as.factor(test.Num),
  Stock=factor(Stock,levels=c("LP","UP"),ordered=F),
  Temp.factor=factor(Temp.factor),
  Pond=factor(all_data$Pond,levels=c("DB","BC","ME"),ordered=F)
)

#############################################
##Functions needed for analysis
#############################################
smooth_diff<-function(model, newdata, f1, f2, var, alpha = 0.05,
                      unconditional = FALSE) {
  xp <- predict(model, newdata = newdata, type = 'lpmatrix')
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  ## zero out cols of X related to splines for other lochs
  X[, ! (c1 | c2)] <- 0
  ## zero out the parametric cols
  X[, !grepl('^s\\(', colnames(xp))] <- 0
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
  crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  data.frame(pair = paste(f1, f2, sep = '-'),
             diff = dif,
             se = se,
             upper = upr,
             lower = lwr)
} ###Taken from Gavin Simpson's blog https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/


smooth_diff.2<-function(model, newdata, f1, f2, d1, var, alpha = 0.05,
                        unconditional = FALSE) {
  xp <- predict(model, newdata = newdata, type = 'lpmatrix')
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  ## zero out cols of X related to splines for other lochs
  X[, ! (c1 | c2)] <- 0
  ## zero out the parametric cols
  X[, !grepl('^s\\(', colnames(xp))] <- 0
  X[, !grepl(d1, colnames(xp))] <- 0
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
  crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  data.frame(pair = paste(f1, f2, sep = '-'),
             diff = dif,
             se = se,
             upper = upr,
             lower = lwr)
} ###Slightly altered version to be able to handle multiple different smoothing parameters needed for MMR




#############################################
###Does the order of testing matter
#############################################
##Is the order of the experiment important for acute repeated measure tests?
#First create a dataframe of repeated temperature measures
Reps<-filter(all_data,all_data$Trial=="Acute"&all_data$test.Num!="2"&all_data$test.Num!="3"&all_data$Experiment<7)
Reps
sapply(Reps,class)
Reps$Experiment<-factor(Reps$Experiment)
Reps$test.Num<-factor(Reps$test.Num)

reps.21<-filter(Reps,Temp.factor==21)
reps.25<-filter(Reps,Temp.factor==25)
#Visualize all of the trials at the same temperature
ggplot(data=reps.21, aes(x=test.Num, y=SMR.abs, color=as.factor(Temp.factor))) + geom_boxplot() + ggtitle("Reps")
ggplot(data=reps.21, aes(x=test.Num, y=MR.Best.abs, color=as.factor(Temp.factor))) + geom_boxplot() + ggtitle("Reps")


ggplot(data=reps.25, aes(x=test.Num, y=SMR.abs, color=as.factor(Temp.factor))) + geom_boxplot() + ggtitle("Reps")
ggplot(data=reps.25, aes(x=test.Num, y=MR.Best.abs, color=as.factor(Temp.factor))) + geom_boxplot() + ggtitle("Reps")

#Test for significant differences among all three for SMR and MMR
reps.test<-aov(SMR.abs~test.Num,data=reps.21)
anova(reps.test)
TukeyHSD(reps.test)

reps.test<-aov(SMR.abs~test.Num,data=reps.25)
anova(reps.test)
TukeyHSD(reps.test)

reps.test<-aov(MR.Best.abs~test.Num,data=reps.21)
anova(reps.test)
TukeyHSD(reps.test)

reps.test<-aov(MR.Best.abs~test.Num,data=reps.25)
anova(reps.test)
TukeyHSD(reps.test)

#SMR is significantly different between Rep 1 and Rep 2 for SMR only lets see if its a specific experiment or not

#Visualize data
ggplot(data=reps.21, aes(x=as.factor(Replicate), y=SMR.abs)) + geom_boxplot() + ggtitle("Acute, SMR 21C")
ggplot(data=reps.21, aes(x = as.factor(Replicate), y = SMR.abs, group=Ind, color=Pond))+geom_line()+geom_point() + ggtitle("Acute, SMR 25C") + theme(legend.position="none")


ggplot(data=reps.25, aes(x=as.factor(Replicate), y=SMR.abs)) + geom_boxplot() + ggtitle("Acute, SMR 21C")
ggplot(data=reps.25, aes(x = as.factor(Replicate), y = SMR.abs, group=Ind, color=Pond))+geom_line()+geom_point() + ggtitle("Acute, SMR 25C") + theme(legend.position="none")


exp.test<-aov(SMR.abs~as.factor(Replicate)*Experiment,data=reps.21)
anova(exp.test)
TukeyHSD(exp.test)

#Experiment 5 is the only experiment with significant differences among rep 1 and rep 2
Exp.5<-filter(reps.21,reps.21$Experiment==5)
ggplot(data=Exp.5, aes(x = as.factor(Replicate), y = SMR.units, group=Ind, color=Ind))+geom_line()+geom_point() + ggtitle("Acute, SMR 21") + theme(legend.position="none")
#mixed bag of individual responses but for fish that are significantly different the trend is to decrease in SMR which is associated with exhaustion.

exp.test<-aov(SMR.abs~as.factor(Replicate)*Experiment,data=reps.25)
anova(exp.test)
TukeyHSD(exp.test)
#No significant differences among experiments


#So only 1 experiment saw a significant different between rep 1 and 2, see manuscript for explanation and rational for only using rep 1
final_data<-filter(all_data,all_data$Replicate==1)

############################################
###Data analysis prep
############################################
final_data$AS.Best.abs<-final_data$MR.Best.abs-final_data$SMR.abs
final_data$AS.Temp.true<-(final_data$SMR.Temp.true+final_data$MMR.Temp.true)/2

#Split data into relevant groupings
acute_spring<-filter(final_data,Trial=="Acute")
acclim_spring<-filter(final_data,Trial=="Acclimated")


###########################################
###General Relationships
##########################################

#Relationship of Three metabolic traits with Mass
ggplot(data=final_data, aes(log(Mass),y=log(SMR.abs),color=Pond))+geom_point(size=2)+labs(x="Mass",y="SMR")+geom_smooth(method="lm")+theme_classic()+theme(legend.position = "None",axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))

ggplot(data=final_data, aes(log(Mass),y=log(MR.Best.abs),color=Pond))+geom_point(size=2)+labs(x="Mass",y="SMR")+geom_smooth(method="lm")+theme_classic()+theme(legend.position = "None",axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))


#Considerations for dealing with the allometric relationship between mass and metabolic rates
# 1. In the acute trials we have repeated measures so we should only look at one temperature at time
# 2. Temperature also has an effect on metabolic rates so allometric scaling should only be considered for one temperature at a time.
# 3. Ponds particularly if they are genetically distinct may have their own scaling as well
# 4. As fish grow and enter different life stages the scaling factor can also change 

all_21<-filter(final_data,Temp.factor=="21")
all_21_acute<-filter(final_data,Temp.factor=="21"&Trial=="Acute")
all_21_acclimated<-filter(final_data,Temp.factor=="21"&Trial=="Acclimated")

all_23<-filter(final_data,Temp.factor=="23")
all_23_acute<-filter(final_data,Temp.factor=="23"&Trial=="Acute")
all_23_acclimated<-filter(final_data,Temp.factor=="23"&Trial=="Acclimated")

all_25<-filter(final_data,Temp.factor=="25")
all_25_acute<-filter(final_data,Temp.factor=="25"&Trial=="Acute")
all_25_acclimated<-filter(final_data,Temp.factor=="25"&Trial=="Acclimated")


mass.smr.21<-lm(log(SMR.abs)~log(Mass),data=all_21_acclimated)
summary(mass.smr.21)
mass.smr.23<-lm(log(SMR.abs)~log(Mass),data=all_23_acclimated)
summary(mass.smr.23)
mass.smr.25<-lm(log(SMR.abs)~log(Mass),data=all_25_acclimated)
summary(mass.smr.25)

#Scaling varies dramatically across the temperatures the 21C scaling factor is the closest to other values seen in the literature .8-.9
#Now lets consider acute vs acclimated and ponds looking at 21C data

mass.trial<-lm(log(SMR.abs)~log(Mass)*Trial,data=all_21)
summary(mass.trial)
mass.pond<-lm(log(SMR.abs)~log(Mass)*Pond,data=all_21)
summary(mass.pond)

#No difference among trials but ME is significantly different than the others
ggplot(data=all_21, aes(log(Mass),y=log(SMR.abs),color=Pond))+geom_point(size=2)+labs(x="Mass",y="SMR")+geom_smooth(method="lm")+theme_classic()+theme(legend.position = "None",axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
ggplot(data=all_21, aes(log(Mass),y=log(MR.Best.abs),color=Pond))+geom_point(size=2)+labs(x="Mass",y="SMR")+geom_smooth(method="lm")+theme_classic()+theme(legend.position = "None",axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))

#So if we want to scale metabolic rates to a similar weight fish as is commonly done in the literature ME fish will have a different scaling factor than DB and BC fish.

ME_21<-filter(all_21,Pond=="ME")#Menominee fish
DB_21<-filter(all_21,Pond=="DB")
BC_21<-filter(all_21,Pond=="BC")

mass.me<-lm(log(SMR.abs)~log(Mass),data=ME_21)
summary(mass.me)
mass.bc<-lm(log(SMR.abs)~log(Mass),data=BC_21)
summary(mass.bc)
mass.db<-lm(log(SMR.abs)~log(Mass),data=DB_21)
summary(mass.db)

mean(final_data$Mass)

for (i in 1:nrow(final_data)){
  if(final_data$Pond[i]=="DB"){
    final_data$SMR.ma[i]<-final_data$SMR.units[i]*(final_data$Mass[i]/10)**(1-0.71)
  } else if (final_data$Pond[i]=="BC") {
    final_data$SMR.ma[i]<-final_data$SMR.units[i]*(final_data$Mass[i]/10)**(1-0.74)
  } else if (final_data$Pond[i]=="ME") {
    final_data$SMR.ma[i]<-final_data$SMR.units[i]*(final_data$Mass[i]/10)**(1-0.82)
  }}

for (i in 1:nrow(all_21)){
  if(all_21$Pond[i]=="DB"){
    all_21$SMR.ma[i]<-all_21$SMR.units[i]*(all_21$Mass[i]/10)**(1-.71)
  } else if (all_21$Pond[i]=="BC") {
    all_21$SMR.ma[i]<-all_21$SMR.units[i]*(all_21$Mass[i]/10)**(1-.74)
  } else if (all_21$Pond[i]=="ME") {
    all_21$SMR.ma[i]<-all_21$SMR.units[i]*(all_21$Mass[i]/10)**(1-.82)
  }}
ggplot(data=all_21, aes(Mass,y=SMR.ma,color=Pond))+geom_point(size=2)+labs(x="Mass",y="SMR")+geom_smooth(method="lm")+theme_classic()+theme(legend.position = ,axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))



#MMR
mass.mmr.21<-lm(log(MR.Best.abs)~log(Mass),data=all_21)
summary(mass.mmr.21)
mass.mmr.23<-lm(log(MR.Best.abs)~log(Mass),data=all_23)
summary(mass.mmr.23)
mass.mmr.25<-lm(log(MR.Best.abs)~log(Mass),data=all_25)
summary(mass.mmr.25)

#Scaling varies dramatically across the temperatures the 21C scaling factor is the closest to other values seen in the literature .8-.9
#Now lets consider acute vs acclimated and ponds looking at 21C data


mass.pond<-lm(log(MR.Best.abs)~log(Mass)*Pond,data=all_21)
summary(mass.pond)
mass.me<-lm(log(MR.Best.abs)~log(Mass),data=ME_21)
summary(mass.me)
mass.bc<-lm(log(MR.Best.abs)~log(Mass),data=BC_21)
summary(mass.bc)
mass.db<-lm(log(MR.Best.abs)~log(Mass),data=DB_21)
summary(mass.db)


for (i in 1:nrow(final_data)){
  if(final_data$Pond[i]=="DB"){
    final_data$MMR.Best.ma[i]<-final_data$MMR.Best.units[i]*(final_data$Mass[i]/10)**(1-0.77143)
  } else if (final_data$Pond[i]=="BC") {
    final_data$MMR.Best.ma[i]<-final_data$MMR.Best.units[i]*(final_data$Mass[i]/10)**(1-0.86142)
  } else if(final_data$Pond[i]=="ME") {
    final_data$MMR.Best.ma[i]<-final_data$MMR.Best.units[i]*(final_data$Mass[i]/10)**(1-0.85873)
  }}

all_21.final<-filter(final_data,Temp.factor=="21")
ggplot(data=final_data, aes(log(Mass),y=log(MMR.Best.ma),color=Pond))+geom_point(size=2)+labs(x="Mass",y="SMR")+geom_smooth(method="lm")+theme_classic()+theme(legend.position = "None",axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))

mass.final<-lm(log(SMR.ma)~log(Mass)*Pond,data=all_21.final)
summary(mass.final)

mmr.mass.final<-lm(log(MMR.Best.ma)~log(Mass)*Pond,data=all_21.final)
summary(mmr.mass.final)

#AS
mass.as.21<-lm(log(AS.Best.abs)~log(Mass)*Pond,data=all_21)
summary(mass.as.21)
ggplot(data=all_21.final, aes(log(Mass),y=log((AS.Best.absolute.units*(Mass/10)**(1-0.88820))),color=Pond))+geom_point(size=2)+labs(x="Mass",y="SMR")+geom_smooth(method="lm")+theme_classic()+theme(legend.position = "None",axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))

mass.me<-lm(log(AS.Best.abs)~log(Mass),data=ME_21)
summary(mass.me)
mass.bc<-lm(log(AS.Best.abs)~log(Mass),data=BC_21)
summary(mass.bc)
mass.db<-lm(log(AS.Best.abs)~log(Mass),data=DB_21)
summary(mass.db)

for (i in 1:nrow(final_data)){
  if(final_data$Pond[i]=="DB"){
    final_data$AS.ma[i]<-final_data$AS.Best.absolute.units[i]*(final_data$Mass[i]/10)**(1-0.80601)
  } else if (final_data$Pond[i]=="BC") {
    final_data$AS.ma[i]<-final_data$AS.Best.absolute.units[i]*(final_data$Mass[i]/10)**(1-0.93737)
  } else if(final_data$Pond[i]=="ME") {
    final_data$AS.ma[i]<-final_data$AS.Best.absolute.units[i]*(final_data$Mass[i]/10)**(1-0.88265)
  }}


final_data$AS.ma<-final_data$MMR.Best.ma-final_data$SMR.ma
ggplot(data=final_data, aes(log(Mass),y=log(AS.ma)))+geom_point(size=2)+labs(x="Mass",y="AS")+geom_smooth(method="lm")+theme_classic()+theme(legend.position = "None",axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))

final_data$AS.fact<-(final_data$MR.Best.abs)/final_data$SMR.abs

#Relationship of SMR with MMR

smrvmmr.lm<-lme(log(MMR.abs)~log(Mass)+MMR.Temp.true*Trial,random=~1|Ind,data=final_data)
summary(smrvmmr.lm)
res.mass<-smrvmmr.lm$residuals[,1]
final_data.2<-final_data
final_data.2$res.mass<-res.mass
ggplot(final_data.2,aes(x=res.mass,y=log(SMR.abs)))+geom_smooth(method=lm)+geom_point()
mass.lm<-lm(res.mass~log(SMR.abs),data=final_data.2)
summary(mass.lm)



smrvmmr.lm<-lme(log(MMR.abs)~log(Mass)+MMR.Temp.true,random=~1|Ind,data=acute_spring)
summary(smrvmmr.lm)
res.mass<-smrvmmr.lm$residuals[,1]
acute_spring.2<-acute_spring
acute_spring.2$res.mass<-res.mass
acute_spring.2<-filter(acute_spring.2,Temp.factor!="27")
acute_spring.2<-filter(acute_spring.2,Temp.factor!="29")
ggplot(acute_spring.2,aes(x=res.mass,y=log(SMR.abs)))+geom_smooth(method=lm)+geom_point()
mass.lm<-lm(res.mass~log(SMR.abs),data=acute_spring.2)
summary(mass.lm)


smrvmmr.lm<-lme(log(MMR.abs)~log(Mass)+MMR.Temp.true,random=~1|Ind,data=acclim_spring)
summary(smrvmmr.lm)
res.mass<-smrvmmr.lm$residuals[,1]
acclim_spring.2<-acclim_spring
acclim_spring.2$res.mass<-res.mass
ggplot(acclim_spring.2,aes(x=res.mass,y=log(SMR.abs)))+geom_smooth(method=lm)+geom_point()
mass.lm<-lm(res.mass~log(SMR.abs),data=acclim_spring.2)
summary(mass.lm)

smrvmmr.pond<-ggplot(final_data, aes(x = log(SMR.ma), y = log(MMR.Best.ma), color=Pond)) + 
  theme_classic()+ geom_point() +
  theme(legend.box.background = element_rect(colour = "black"),axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),plot.title = element_text(size=25,face="bold",hjust = 0.5))+
  guides(color=guide_legend(title="Pond"),shape=guide_legend(title="Population")) +
  scale_shape_manual(values=c(16,17,15),labels=c('Dearborn','Bay City','Menominee'))+geom_smooth(method=lm)
smrvmmr.pond


smrvmmr.temp<-ggplot(final_data, aes(x = log(SMR.ma), y = log(MMR.Best.ma), color=Temp.factor)) + 
  theme_classic()+ geom_point() +
  theme(legend.box.background = element_rect(colour = "black"),axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),plot.title = element_text(size=25,face="bold",hjust = 0.5))+
  guides(color=guide_legend(title="Pond"),shape=guide_legend(title="Population")) +
  scale_shape_manual(values=c(16,17,15),labels=c('Dearborn','Bay City','Menominee'))+geom_smooth(method=lm)
smrvmmr.temp

smrvmmr.mass<-ggplot(final_data, aes(x = SMR.ma, y = MR.Best.ma, color=Mass)) + 
  theme_classic()+ geom_point() + scale_y_continuous(expand = expansion(mult = c(0, 0)),limits=c(0,20)) + 
  labs(y="",x=NULL) 

#theme(legend.box.background = element_rect(colour = "black"),axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),plot.title = element_text(size=25,face="bold",hjust = 0.5))+
#guides(color=guide_legend(title="Temperature"),shape=guide_legend(title="Population")) + scale_color_manual(values=c("blue","lightblue","gold","darkgoldenrod1","darkred"))+
#scale_shape_manual(values=c(16,17,15),labels=c('Dearborn','Bay City','Menominee'))+geom_smooth(method="gam")



grid.arrange(
  smrvmmr.pond,
  smrvmmr.temp,
  smrvmmr.mass,
  ncol=3,
  top=textGrob("SMR vs MMR",gp=gpar(fontsize=20,font=2)),
  bottom=textGrob("SMR (mgO2/min)",gp=gpar(fontsize=20,font=2)),
  padding=unit(2,"line"))




#Split data into relevant groupings

write.csv(final_data,"final_data.csv")
final_data<-read.csv("final_data.csv")
final_data<- transform(
  final_data,
  Ind=as.factor(Ind),
  SMR.units=as.numeric(SMR.units),
  test.Num=as.factor(test.Num),
  Stock=factor(Stock,levels=c("LP","UP"),ordered=F),
  Temp.factor=factor(Temp.factor),
  Pond=factor(final_data$Pond,levels=c("DB","BC","ME"),ordered=FALSE)
)

acute_spring<-filter(final_data,Trial=="Acute")
acute_spring_all<-filter(all_data,Trial=="Acute")
acclim_spring<-filter(final_data,Trial=="Acclimated")


#########################################
###Data Summary Information
#########################################
DB.acute<-filter(acute_spring,Pond=="DB")
BC.acute<-filter(acute_spring,Pond=="BC")
ME.acute<-filter(acute_spring,Pond=="ME")
min(DB.acute$Mass)
max(DB.acute$Mass)
mean(DB.acute$Mass)
min(BC.acute$Mass)
max(BC.acute$Mass)
mean(BC.acute$Mass)
min(ME.acute$Mass)
max(ME.acute$Mass)
mean(ME.acute$Mass)

DB.acclim.21<-filter(acclim_spring,Pond=="DB" &Temp.factor=="21")
BC.acclim.21<-filter(acclim_spring,Pond=="BC" &Temp.factor=="21")
ME.acclim.21<-filter(acclim_spring,Pond=="ME" &Temp.factor=="21")
DB.acclim.23<-filter(acclim_spring,Pond=="DB" &Temp.factor=="23")
BC.acclim.23<-filter(acclim_spring,Pond=="BC" &Temp.factor=="23")
ME.acclim.23<-filter(acclim_spring,Pond=="ME" &Temp.factor=="23")
DB.acclim.25<-filter(acclim_spring,Pond=="DB" &Temp.factor=="25")
BC.acclim.25<-filter(acclim_spring,Pond=="BC" &Temp.factor=="25")
ME.acclim.25<-filter(acclim_spring,Pond=="ME" &Temp.factor=="25")

min(DB.acclim.21$Mass)
max(DB.acclim.21$Mass)
mean(DB.acclim.21$Mass)
min(BC.acclim.21$Mass)
max(BC.acclim.21$Mass)
mean(BC.acclim.21$Mass)
min(ME.acclim.21$Mass)
max(ME.acclim.21$Mass)
mean(ME.acclim.21$Mass)

min(DB.acclim.23$Mass)
max(DB.acclim.23$Mass)
mean(DB.acclim.23$Mass)
min(BC.acclim.23$Mass)
max(BC.acclim.23$Mass)
mean(BC.acclim.23$Mass)
min(ME.acclim.23$Mass)
max(ME.acclim.23$Mass)
mean(ME.acclim.23$Mass)

min(DB.acclim.25$Mass)
max(DB.acclim.25$Mass)
mean(DB.acclim.25$Mass)
min(BC.acclim.25$Mass)
max(BC.acclim.25$Mass)
mean(BC.acclim.25$Mass)
min(ME.acclim.25$Mass)
max(ME.acclim.25$Mass)
mean(ME.acclim.25$Mass)


#########################################
###Standard Metabolic Rate
#########################################

##Investigate non-linear relationships


#SMR_gam.1<-gam(SMR.ma~Pond+s(Mass,by=Pond)+s(SMR.Temp.true,by=Pond,k=5)+s(Ind,bs="re"), data=acute_spring, method="REML",select=T)

SMR_gam.1<-gam(SMR.abs~Pond+s(Mass,k=5)+s(SMR.Temp.true,by=Pond,k=5)+s(Ind,bs="re"), data=acute_spring, method="REML",select=T)


acute_spring2<-filter(acute_spring,Pond!="ME")
acute_spring2<- transform(acute_spring, Pond_in_Stock = interaction(Pond, Stock, drop = TRUE))

draw(SMR_gam.1,residuals=TRUE)
draw(mmr.mod_gam.1,residuals=TRUE)
draw(as.mod_gam.1,residualts=TRUE)
#Rational for structure here: By variable allows the non-linear relationship to vary between the levels of populations,
#However, by factor smooths are centered on 0 so we have to include the parametric term for populations to include the group level means into the model
#the id term ensures that while they can vary they are of similar wiggliness,
#k determines the number of basis functions to be used and here it is set for 5 1 for each temperature factor
#Individual is treated as a random factor
#Method REML 
#Select=T this fully penalizes the smooths allowing gam to reduce functions to 0 if they are not playing any role in the model


SMR_gam.1
gam.check(SMR_gam.1)#Check the k significance for the model
appraise(SMR_gam.1)#Test the model
predict(SMR_gam.1)
plot(SMR_gam.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
visreg(SMR_gam.1)#Visualize relationships
summary(SMR_gam.1)
concurvity(SMR_gam.1) #checking for concurvity amongst smooths, non-parametric form of collinearity test for GAM

concurvity(SMR_gam.1,full=F)
#Seeing some concurvity in the model, Mass, Individual and Pond, however, both of these variables are important due to the allometric relationship of smr and mass as well as the repeated measure

tab_model(SMR_gam.1)

##Test replicates
acute_test<-filter(acute_spring,Experiment!=5)
SMR_gam.rep<-gam(SMR.abs~Pond+s(Mass)+s(SMR.Temp.true,by=Pond,k=5)+s(Ind,bs="re"), data=acute_test, method="REML",select=T)
gam.check(SMR_gam.rep)#Check the k significance for the model
appraise(SMR_gam.rep)#Test the model
plot(SMR_gam.rep, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
visreg(SMR_gam.rep)#Visualize relationships
summary(SMR_gam.rep)
concurvity(SMR_gam.rep) #checking for concurvity amongst smooths, non-parametric form of collinearity test for GAM


#Test to see if outliers from repeated measures are having significant effect
acute_spring2<-filter(acute_spring,acute_spring$Ind!=24)
acute_spring2<-filter(acute_spring,acute_spring$Experiment!=5)
SMR_gam.a<-gam(SMR.abs~Pond+s(Mass)+s(SMR.Temp.true,by=Pond,k=5)+s(Ind,bs="re"), data=acute_spring2, method="REML",select=T)

gam.check(SMR_gam.a)#Check the k significance for the model
appraise(SMR_gam.a)#Test the model
visreg(SMR_gam.a)#Visualize relationships
summary(SMR_gam.a)

#No significant differences in models

#simplified models to assess R2 and deviance
SMR_gam.simple<-gam(SMR.abs~s(Mass), data=acute_spring, method="REML",select=T)
summary(SMR_gam.simple)
appraise(SMR_gam.simple)

SMR_gam.simple.2<-gam(SMR.abs~s(SMR.Temp.true,by=Pond,k=5), data=acute_spring, method="REML",select=T)
summary(SMR_gam.simple.2)
appraise(SMR_gam.simple.2)

#Still have very high R2 and deviance explained even in simplistic models. 

###Test for differences between parametric terms
library(itsadug)

wald_gam(SMR_gam.1)


#Plot of Effect of Populations
smr.pop<-visreg(SMR_gam.1,xvar='Pond',gg=T) + labs(x="",y=expression(paste('Metabolic Rate '(mg*O[2] * min^-1)))) + ggtitle("SMR")+ ylim(0,2.5)+
  theme_classic()+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds"))
smr.pop


#Plot Smooths
smr.gamS<-smooth_estimates(SMR_gam.1)|>
  add_confint()
smr.gamS<-filter(smr.gamS,smr.gamS$.smooth!="s(Mass)")
smr.gamS<-filter(smr.gamS,smr.gamS$.type=="TPRS")
smr.gamS$smooth<-factor(smr.gamS$.smooth,levels=c("s(SMR.Temp.true):PondDB","s(SMR.Temp.true):PondBC","s(SMR.Temp.true):PondME"))
smr.gamS<-smr.gamS %>%
  add_constant(coef(SMR_gam.1)["(Intercept)"]) %>%
  transform_fun(inv_link(SMR_gam.1)) 

coef(SMR_gam.1)
DB.r<-filter(smr.gamS,Pond=="DB")


BC.r<-filter(smr.gamS,Pond=="BC")
BC.r<-BC.r%>%
  add_constant(-0.05)

ME.r<-filter(smr.gamS,Pond=="ME")
ME.r<-ME.r%>%
  add_constant(.229)

all.smr<-rbind(DB.r,BC.r,ME.r)

smr.graph<-ggplot(data=all.smr, aes(x=SMR.Temp.true,y=.estimate,color=Pond)) + 
  geom_line(size=1) + 
  labs(x="",y="") +
  scale_color_manual(values=c("deeppink2","orange","deepskyblue"),labels=c('Dearborn:LP','Bay City:LP','Menominee:UP')) + 
  theme_classic()+theme(legend.position="",legend.text=element_text(size=rel(1)),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds"))+scale_x_continuous(breaks=c(21,23,25,27,29))
smr.graph<-smr.graph+geom_ribbon(aes(ymin=.lower_ci,ymax=.upper_ci,x=SMR.Temp.true,fill=Pond),linetype=0,alpha=.25)+scale_fill_manual(values=c("#FFC0CB80", "#FFD70080","#ADD8E680"))
smr.graph




#Mass smooth
smr.gamS<-smooth_estimates(SMR_gam.1)|>
  add_confint()
smr.gamS<-filter(smr.gamS,smr.gamS$.smooth=="s(Mass)")
smr.gamS<-filter(smr.gamS,smr.gamS$.type=="TPRS")
smr.gamS<-smr.gamS %>%
  add_constant(coef(SMR_gam.1)["(Intercept)"]) %>%
  transform_fun(inv_link(SMR_gam.1)) 


mass.smr.graph<-ggplot(data=all.smr, aes(x=Mass,y=.estimate,color=Pond)) +
  geom_line(size=1) + 
  labs(x="",y=expression(paste('Metabolic Rate '(mg*O[2] * min^-1)))) + ggtitle("SMR")+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue"),labels=c('Dearborn:LP','Bay City:LP','Menominee:UP')) + 
  theme_classic()+theme(legend.position="",legend.text=element_text(size=rel(1)),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds"))
mass.smr.graph<-mass.smr.graph+geom_ribbon(aes(ymin=smr.gamS$.lower_ci,ymax=smr.gamS$.upper_ci,x=smr.gamS$Mass),linetype=0,alpha=.1)
mass.smr.graph

##Check difference between smooths

#First create matrix of values to be plugged into models
pdat<-expand.grid(SMR.Temp.true=seq(21,29,length=400),Mass=mean(acute_spring$Mass),
                  Pond=c('DB',"BC","ME"),Ind=1)
Ind<-1
Ind<-c(Ind,Ind,Ind)
pdat<-cbind(pdat,Ind)

#Using smooth_difff equation do a pairwise comparions among the populations
comp1<-smooth_diff(SMR_gam.1,pdat,'DB','BC','Pond')
comp2<-smooth_diff(SMR_gam.1,pdat,'DB','ME','Pond')
comp3<-smooth_diff(SMR_gam.1,pdat,'BC','ME','Pond')
comp.2<-cbind(SMR.Temp.true=seq(21,29,length=400),
              rbind(comp1,comp2,comp3))
smr.c1<-ggplot(comp.2, aes(x = SMR.Temp.true, y = diff, group = pair)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() + ggtitle("SMR")+
  facet_wrap(~ pair, ncol = 3,labeller="label_both") +
  coord_cartesian(ylim = c(-1,1)) +
  labs(x = "", y = '')+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
smr.c1<-smr.c1+geom_hline(yintercept=0)
smr.c1


smr.pred<-predict(SMR_gam.1,pdat,se.fit=TRUE)
smr.pred<-cbind(smr.pred,pdat)
ilink<-family(SMR_gam.1)$linkinv
smr.pred<-transform(smr.pred,
          lwr_ci=ilink(fit-(2*se.fit)),
          upr_ci=ilink(fit+(2*se.fit)),
          fitted=ilink(fit))
        
ggplot(smr.pred, aes(x = SMR.Temp.true, y = fitted,color=Pond)) +
  geom_ribbon(aes(ymin = lwr_ci, ymax = upr_ci), alpha = 0.1,linetype=0) +
  geom_line()
smr_dif<-difference_smooths(SMR_gam.1,select="s(SMR.Temp.true)",group_means = TRUE)
dif_smr<-draw(smr_dif,ncol=3)
dif_smr<-dif_smr & xlab("Temperature (\u00B0C)") & ylab(expression(paste('Metabolic Rate '(mg*O[2] * min^-1))))
mmr_dif<-difference_smooths(mmr.mod_gam.1,select="s(MMR.Temp.true)",group_means = TRUE)
dif_mmr<-draw(mmr_dif,ncol=3)
dif_mmr<-dif_mmr& xlab("Temperature (\u00B0C)") & ylab(expression(paste('Metabolic Rate '(mg*O[2] * min^-1))))
as_dif<-difference_smooths(as.mod_gam.1,select="s(AS.Temp.true)",group_means = TRUE)
dif_as<-draw(as_dif,ncol=3)
dif_as<-dif_as & xlab("Temperature (\u00B0C)") & ylab(expression(paste('Metabolic Rate '(mg*O[2] * min^-1))))

library(patchwork)
dif_smr / dif_mmr / dif_as 
dif_mmr
dif_as



smr_diff<-ggemmeans(SMR_gam.1,terms=c("SMR.Temp.true","Pond"))
plot(smr_diff)
emmeans(SMR_gam.1,specs=c("SMR.Temp.true","Pond","Mass"),at=list(SMR.Temp.true=seq(from=21,to=29,by=2), Mass=10))
emmeans(SMR_gam.1,pairwise~Pond|SMR.Temp.true,at=list(SMR.Temp.true=seq(from=21,to=29,by=2), Mass=10))
###################################
###Maximum Metabolic Rate
###################################


#mmr.mod_gam.1<-gam(MR.Best.abs~s(SMR.abs,k=5)+s(Mass)+s(Ind,bs="re")+Pond+s(MMR.Temp.true,by=Pond,k=5),data=acute_spring, method="REML",select=T)
mmr.mod_gam.1<-gam(MR.Best.abs~ s(Mass)+s(Ind,bs="re")+Pond+s(MMR.Temp.true,by=Pond,k=5),
                   data=acute_spring, method="REML",select=T)




mmr.mod_gam.1
gam.check(mmr.mod_gam.1)#Check the k significance for the model
visreg(mmr.mod_gam.1)
plot(mmr.mod_gam.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(mmr.mod_gam.1)
summary(mmr.mod_gam.1)
draw(mmr.mod_gam.1)
concurvity(mmr.mod_gam.1)
#Seeing concurivity here too
concurvity(mmr.mod_gam.1,full=F)
#concurvity seems to stem from Mass, Individual and Pond

mmr.mod_gam.2
gam.check(mmr.mod_gam.2)#Check the k significance for the model
visreg(mmr.mod_gam.2)
plot(mmr.mod_gam.2, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(mmr.mod_gam.2)
summary (mmr.mod_gam.2)
draw(mmr.mod_gam.2)
concurvity(mmr.mod_gam.2)
#Seeing concurivity here too
concurvity(mmr.mod_gam.2,full=F)
#concurvity seems to stem from Mass, Individual and Pond

AIC(mmr.mod_gam.1,mmr.mod_gam.2)

wald_gam(mmr.mod_gam.1)

#Single out effect of population for plot
mmr.pop<-visreg(mmr.mod_gam.1,xvar='Pond',gg=T)+ggtitle("MMR")+labs(x="",y=expression(paste('Metabolic Rate '(mg*O[2] * min^-1))))+ylim(0,5)+
  theme_classic()+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds")) 
mmr.pop

#Plot Temperature smooths from model
s.gamr<-smooth_estimates(mmr.mod_gam.1)|>
  add_confint()
t.gamr<-s.gamr[!is.na(s.gamr$MMR.Temp.true),]
t.gamr<-t.gamr %>%
  add_constant(coef(mmr.mod_gam.1)["(Intercept)"]) %>%
  transform_fun(inv_link(mmr.mod_gam.1)) 
 

coef(mmr.mod_gam.1)


DB.r<-filter(t.gamr,Pond=="DB")


BC.r<-filter(t.gamr,Pond=="BC")
BC.r<-BC.r%>%
  add_constant(.217)

ME.r<-filter(t.gamr,Pond=="ME")
ME.r<-ME.r%>%
  add_constant(.729)

all.mmr<-rbind(DB.r,BC.r,ME.r)

mmr.graph<-ggplot(data=all.mmr, aes(x=MMR.Temp.true,y=.estimate,color=Pond)) +
  geom_line(size=1,aes(linetype=Pond)) + scale_linetype_manual(values=c("solid", "solid","solid"))+guides(linetype ="none" )+
  labs(x="",y="") +
  scale_color_manual(values=c("deeppink2","orange","deepskyblue"),labels=c("Dearborn:LP","Bay City:LP","Menominee:UP","Global")) + 
  theme_classic()+theme(legend.position="",legend.text=element_text(size=rel(1)),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Populations"))+
  scale_x_continuous(breaks=c(21,23,25,27,29))
  #geom_vline(xintercept=ME.r, color="deepskyblue",linetype="dashed")
mmr.graph<-mmr.graph+geom_ribbon(aes(ymin=.lower_ci,ymax=.upper_ci,x=MMR.Temp.true,fill=Pond),linetype=0,alpha=.25)+scale_fill_manual(values=c("#FFC0CB80", "#FFD70080","#ADD8E680"))
mmr.graph

#Plot SMR smooths from model
#s.gamr<-smooth_estimates(mmr.mod_gam.1)
#t.gamr<-s.gamr[!is.na(s.gamr$SMR.abs),]
#mmr.st.1<-ggplot(data=t.gamr, aes(x=SMR.abs,y=.estimate)) +
#  geom_line(size=2) +labs(x="SMR (mgO2/kg/min/\u00B0C)",y="SMR Effect Trend") +
#  ggtitle("SMR Effect on MMR")+
#  scale_color_manual(values=c("blue","lightblue","gold","darkgoldenrod1","darkred"),labels=c("21\u00B0C","23\u00B0C","25\u00B0C","27\u00B0C","29\u00B0C","Global")) +
#  theme_classic()+theme(legend.text=element_text(size=rel(1.5)),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Temperature"))
#mmr.st.1


##Mass smooth from MMR
mmr.gamm<-smooth_estimates(mmr.mod_gam.1)|>
  add_confint()
mmr.gamm<-filter(mmr.gamm,mmr.gamm$.smooth=="s(Mass)")
mmr.gamm<-filter(mmr.gamm,mmr.gamm$.type=="TPRS")
mmr.gamm<-mmr.gamm %>%
  add_constant(coef(mmr.mod_gam.1)["(Intercept)"]) %>%
  transform_fun(inv_link(SMR_gam.1)) 
mass.mmr.graph<-ggplot(data=mmr.gamm, aes(x=Mass,y=.estimate,color=Pond)) +
  geom_line(size=2) + 
  labs(x="Temperature (\u00B0C)",y=expression(paste('Metabolic Rate '(mg*O[2] * min^-1)))) + ggtitle("MMR")+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue"),labels=c('Dearborn:LP','Bay City:LP','Menominee:UP')) + 
  theme_classic()+theme(legend.position="",legend.text=element_text(size=rel(1)),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds"))
mass.mmr.graph<-mass.mmr.graph+geom_ribbon(aes(ymin=mmr.gamm$.lower_ci,ymax=mmr.gamm$.upper_ci,x=mmr.gamm$Mass),linetype=0,alpha=.1)
mass.mmr.graph

#Check difference between smooths

#First create matrix of values to be plugged into models
pdat<-expand.grid(MMR.Temp.true=seq(21,29,length=400),Mass=seq(min(acute_spring$Mass),max(acute_spring$Mass),length=400),
                  Pond=c('DB',"BC","ME"))

SMR.abs<-seq(min(acute_spring$SMR.abs),max(acute_spring$SMR.abs),length=400)
SMR.abs<-c(SMR.abs,SMR.abs,SMR.abs)
Ind<-seq(1,400)
Ind<-c(Ind,Ind,Ind)
temp1<-rep(21, length=80)
temp2<-rep(23, length=80)
temp3<-rep(25, length=80)
temp4<-rep(27, length=80)
temp5<-rep(29, length=80)
Temp.factor<-c(temp1,temp2,temp3,temp4,temp5)


pdat<-cbind(pdat,Ind,SMR.abs,Temp.factor)


#Using smooth_diff.2 which allows us to specify which smooths we are interested in comparing
#First temperature smooths
comp1<-smooth_diff.2(mmr.mod_gam.1,pdat,'DB','BC','MMR.Temp.true','Pond')
comp2<-smooth_diff.2(mmr.mod_gam.1,pdat,'DB','ME','MMR.Temp.true','Pond')
comp3<-smooth_diff.2(mmr.mod_gam.1,pdat,'BC','ME','MMR.Temp.true','Pond')
comp.mmr.1<-cbind(MMR.Temp.true=seq(21,29,length=400),
                  rbind(comp1,comp2,comp3))
head(comp.mmr.1)
mmr.c.1<-ggplot(comp.mmr.1, aes(x = MMR.Temp.true, y = diff, group = pair)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +ggtitle("MMR")+
  facet_wrap(~ pair, ncol = 3,labeller="label_both") +
  coord_cartesian(ylim = c(-1,1)) +
  labs(x = "", y = 'Difference in Temperature Trend')+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
mmr.c.1<-mmr.c.1+geom_hline(yintercept=0)
mmr.c.1

#######################################
###Aerobic Scope
#######################################

as.mod_gam.1<-gam(AS.Best.abs~Pond+s(Mass,k=5)+s(AS.Temp.true,by=Pond,k=5)+#allow effect of temp factor to have varying intercepts
                     s(Ind,bs="re"),#random effect for individual fish
                   data=acute_spring, method="REML",select=T)


gam.check(as.mod_gam.1)
summary(as.mod_gam.1)
plot(as.mod_gam.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
anova(as.mod_gam.1)
draw(as.mod_gam.1)
appraise(as.mod_gam.1)
visreg(as.mod_gam.1)
draw(as.mod_gam.1)

#Graph Population Effect

wald_gam(as.mod_gam.1)

as.pop<-visreg(as.mod_gam.1,xvar='Pond',gg=T)+labs(x="Pond Populations",y=expression(paste('Metabolic Rate '(mg*O[2] * min^-1)))) + ggtitle("AS")+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue")) + ylim(0,3.5) +
  theme_classic()+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds")) 
as.pop

coef(as.mod_gam.1)["(Intercept)"]

#Graph smooths
s.gama<-smooth_estimates(as.mod_gam.1)|>
  add_confint(coverage=0.95)
s.gama<-s.gama %>%
  add_constant(coef(as.mod_gam.1)["(Intercept)"]) %>%
  transform_fun(inv_link(as.mod_gam.1)) 
s.gama<-filter(s.gama,s.gama$.type=="TPRS")
s.gama<-s.gama[!is.na(s.gama$AS.Temp.true),]
s.gama$smooth<-factor(s.gama$.smooth,levels=c("s(AS.Temp.true):PondDB","s(AS.Temp.true):PondBC","s(AS.Temp.true):PondME"))
as.acute<- acute_spring|>
  add_partial_residuals(as.mod_gam.1)

coef(as.mod_gam.1)

DB.r<-filter(s.gama,Pond=="DB")



DB.max<-max(DB.r$.estimate)
DB.max<-filter(DB.r,DB.r$.estimate==DB.max)
DB.max<-DB.max$AS.Temp.true

BC.r<-filter(s.gama,Pond=="BC")
BC.r<-BC.r%>%
  add_constant(0.247)



BC.max<-max(BC.r$.estimate)
BC.max<-filter(BC.r,BC.r$.estimate==BC.max)
BC.max<-BC.max$AS.Temp.true

ME.r<-filter(s.gama,Pond=="ME")
ME.r<-ME.r%>%
  add_constant(.522)

ME.max<-max(ME.r$.estimate)
ME.max<-filter(ME.r,ME.r$.estimate==ME.max)
ME.max<-ME.max$AS.Temp.true

all.as<-rbind(DB.r,BC.r,ME.r)


as.graph<-ggplot(data=all.as, aes(x=AS.Temp.true,y=.estimate,color=smooth)) +
  geom_line(size=1,aes(linetype=Pond)) + scale_linetype_manual(values=c("solid", "twodash","solid"))+guides(linetype ="none" )+
  labs(x="Temperature (\u00B0C)",y="") +scale_x_continuous(breaks=c(21,23,25,27,29))+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue"),labels=c("Dearborn:LP","Bay City:LP","Menominee:UP")) + 
  theme_classic()+theme(legend.position="",legend.text=element_text(size=rel(1)),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Populations"))#+geom_vline(xintercept=ME.r, color="deepskyblue",linetype="dashed")+geom_vline(xintercept=BC.r, color="orange",linetype="dashed")
as.graph<-as.graph+geom_ribbon(aes(ymin=.lower_ci,ymax=.upper_ci,x=AS.Temp.true,fill=Pond),linetype=0,alpha=.25)+scale_fill_manual(values=c("#FFC0CB80", "#FFD70080","#ADD8E680"))
as.graph





#AS mass smooth
as.gamm<-smooth_estimates(as.mod_gam.1)|>
  add_confint()
as.gamm<-filter(as.gamm,as.gamm$.smooth=="s(Mass)")
as.gamm<-filter(as.gamm,as.gamm$.type=="TPRS")
as.gamm<-as.gamm %>%
  add_constant(coef(as.mod_gam.1)["(Intercept)"]) %>%
  transform_fun(inv_link(SMR_gam.1)) 
mass.as.graph<-ggplot(data=as.gamm, aes(x=Mass,y=.estimate,color=Pond)) +
  geom_line(size=1) + 
  labs(x="",y=expression(paste('Metabolic Rate '(mg*O[2] * min^-1)))) + ggtitle("AS")+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue"),labels=c('Dearborn:LP','Bay City:LP','Menominee:UP')) + 
  theme_classic()+theme(legend.position="",legend.text=element_text(size=rel(1)),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds"))
mass.as.graph<-mass.as.graph+geom_ribbon(aes(ymin=as.gamm$.lower_ci,ymax=as.gamm$.upper_ci,x=as.gamm$Mass),linetype=0,alpha=.1)
mass.as.graph






#Check differences among smooths
pdat<-expand.grid(AS.Temp.true=seq(21,29,length=400),Mass=(10),Pond=c('BC','DB','ME'))
Ind<-seq(1,400)
Ind<-c(Ind,Ind,Ind)
pdat<-cbind(pdat,Ind)
comp1 <- smooth_diff(as.mod_gam.1, pdat, 'DB', 'BC', 'Pond')
comp2 <- smooth_diff(as.mod_gam.1, pdat, 'DB', 'ME', 'Pond')
comp3 <- smooth_diff(as.mod_gam.1, pdat, 'BC', 'ME', 'Pond')
comp.as.2 <- cbind(AS.Temp.true=seq(21,29,length=400),Pond=c('BC','DB','ME'),rbind(comp1, comp2, comp3))

as.comp<-ggplot(comp.as.2, aes(x = AS.Temp.true, y = diff, group = pair)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() + ggtitle("AS")+
  facet_wrap(~ pair, ncol = 3,labeller="label_both") +
  coord_cartesian(ylim = c(-1,1)) + labs(x = "Temperature (\u00B0C)", y = '')+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))

as.comp<-as.comp+geom_hline(yintercept=0)
as.comp





#################################
#Other Graphs
#################################
grid.arrange(
  smr.pop,
  smr.graph,
  mmr.pop,
  mmr.graph,
  as.pop,
  as.graph,
  ncol=2
)

grid.arrange(
  smr.c1+scale_x_continuous(breaks=seq(21,29,by=2)),
  mmr.c.1+scale_x_continuous(breaks=seq(21,29,by=2)),
  as.comp+scale_x_continuous(breaks=seq(21,29,by=2)),
  ncol=1
)

grid.arrange(
  mass.smr.graph,
  mass.mmr.graph,
  mass.as.graph,
  ncol=3
)


#Acute Metabolic Phenotypes
smr.preds<-predict(SMR_gam.1,se.fit=TRUE)
smr_data<-data.frame(smr=smr.preds$fit,low=(smr.preds$fit-2*smr.preds$se.fit),high=(smr.preds$fit+2*smr.preds$se.fit))
smr_data$temp<-acute_spring$SMR.Temp.true
smr_data$SMR.abs<-acute_spring$SMR.abs
smr_data$Pond<-acute_spring$Pond
smr.graph.full<-ggplot(smr_data,aes(x=temp))+
  geom_smooth(aes(y=smr ,color=Pond), size=1,method="gam",se=FALSE, alpha=.5,stat="smooth")+
  geom_point(aes(y=SMR.abs,color=Pond))+theme_classic()+ylab(expression(paste(bold('Metabolic Rate '(mg*O[2] * min^-1)))))+xlab(expression(paste(bold("Temperature (\u00B0C)"))))+scale_color_manual(values=c("deeppink2","orange","deepskyblue"),name="Rearing Pond",labels=c('Dearborn:LP','Bay City:LP','Menominee:UP'))+
  ggtitle(expression(paste(bold("Standard Metabolic Rate"))))+theme(plot.title = element_text(hjust=0.5),legend.position="none")
smr.graph.full<-smr.graph.full+geom_ribbon(aes(ymin=low,ymax=high,y=smr,fill=Pond),stat="smooth",alpha=0.2)+scale_fill_manual(values=c("#FFC0CB80", "#FFD70080","#ADD8E680"),name="Rearing Pond",labels=c('Dearborn:LP','Bay City:LP','Menominee:UP'))+scale_x_continuous(breaks=c(21,23,25,27,29))
smr.graph.full

mmr.preds<-predict(mmr.mod_gam.1,se.fit=TRUE)
mmr_data<-data.frame(mmr=mmr.preds$fit,low=(mmr.preds$fit-2*mmr.preds$se.fit),high=(mmr.preds$fit+2*mmr.preds$se.fit))
mmr_data$temp<-acute_spring$MMR.Temp.true
mmr_data$MR.Best.abs<-acute_spring$MR.Best.abs
mmr_data$Pond<-acute_spring$Pond
mmr.graph.full<-ggplot(mmr_data,aes(x=temp))+
  geom_smooth(aes(y=mmr ,color=Pond), size=1,method="gam",se=FALSE, alpha=.5,stat="smooth")+
  geom_point(aes(y=MR.Best.abs,color=Pond))+theme_classic()+ylab(expression(paste(bold('Metabolic Rate '(mg*O[2] * min^-1)))))+xlab(expression(paste(bold("Temperature (\u00B0C)"))))+scale_color_manual(values=c("deeppink2","orange","deepskyblue"),name="Rearing Pond",labels=c('Dearborn:LP','Bay City:LP','Menominee:UP'))+
  ggtitle(expression(paste(bold("Maximum Metabolic Rate"))))+theme(plot.title = element_text(hjust=0.5),legend.position="none")
mmr.graph.full<-mmr.graph.full+geom_ribbon(aes(ymin=low,ymax=high,y=mmr,fill=Pond),stat="smooth",alpha=0.2)+scale_fill_manual(values=c("#FFC0CB80", "#FFD70080","#ADD8E680"),name="Rearing Pond",labels=c('Dearborn:LP','Bay City:LP','Menominee:UP'))+scale_x_continuous(breaks=c(21,23,25,27,29))
mmr.graph.full

as.preds<-predict(as.mod_gam.1,se.fit=TRUE)
as_data<-data.frame(as=as.preds$fit,low=(as.preds$fit-2*as.preds$se.fit),high=(as.preds$fit+2*as.preds$se.fit))
as_data$temp<-acute_spring$AS.Temp.true
as_data$AS.Best.abs<-acute_spring$AS.Best.abs
as_data$Pond<-acute_spring$Pond
as.graph.full<-ggplot(as_data,aes(x=temp))+
  geom_smooth(aes(y=as ,color=Pond), size=1,method="gam",se=FALSE, alpha=.5,stat="smooth")+
  geom_point(aes(y=AS.Best.abs,color=Pond))+theme_classic()+ylab(expression(paste(bold('Metabolic Rate '(mg*O[2] * min^-1)))))+xlab(expression(paste(bold("Temperature (\u00B0C)"))))+scale_color_manual(values=c("deeppink2","orange","deepskyblue"),name="Rearing Pond",labels=c('Dearborn:LP','Bay City:LP','Menominee:UP'))+
  ggtitle(expression(paste(bold("Aerobic Scope"))))+theme(plot.title = element_text(hjust=0.5),legend.position="none")
as.graph.full<-as.graph.full+geom_ribbon(aes(ymin=low,ymax=high,y=as,fill=Pond),stat="smooth",alpha=0.3)+scale_fill_manual(values=c("#FFC0CB80", "#FFD70080","#ADD8E680"),name="Rearing Pond",labels=c('Dearborn:LP','Bay City:LP','Menominee:UP'))+ scale_x_continuous(breaks=c(21,23,25,27,29))
as.graph.full

grid.arrange(
  smr.graph.full,
  mmr.graph.full,
  as.graph.full,
  ncol=3
)




##No longer using this graph
#grid.arrange(
#  ggplot(acute_spring, aes(x = as.factor(Temp.factor), y = SMR.ma, color=Pond)) + geom_boxplot(outlier.shape = NA,lwd=.5)  + 
#   theme_classic() + scale_y_continuous(expand = expansion(mult = c(0, 0)),limits=c(0,20)) + geom_point(position = position_jitterdodge())+
#   labs(y=NULL,x=NULL)+scale_color_manual(values=c("deeppink2","orange","deepskyblue"),labels=c('Dearborn, LP','Bay City, LP','Menominee, UP'))+
#   ggtitle("Standard Metabolic Rate")+ theme(legend.position=c(.3,.8),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Populations")),
# ggplot(acute_spring, aes(x = as.factor(Temp.factor), y = MMR.Best.ma, color=Pond)) +
#   geom_boxplot(outlier.shape = NA,lwd=.5)  + theme_classic() + 
#   scale_y_continuous(expand = expansion(mult = c(0, 0)),limits=c(0,20)) + 
#   geom_point(position = position_jitterdodge())+
#   labs(y=NULL,x=NULL)+
#   scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+ggtitle("Maximum Metabolic Rate")+
#   theme(legend.position="none",axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds")),
# ggplot(acute_spring, aes(x = as.factor(Temp.factor), y = AS.ma, color=Pond)) + 
#   geom_boxplot(outlier.shape = NA,lwd=.5)  + theme_classic() + 
#   scale_y_continuous(expand = expansion(mult = c(0, 0)),limits=c(0,20)) + 
#   geom_point(position = position_jitterdodge())+
#   labs(y=NULL,x=NULL)+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+ggtitle("Aerobic Scope")+
#   theme(legend.position="none",axis.text=element_text(size=16),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds")),
# ncol=3,
# left= textGrob("Metabolic Rate (mgO2/kg/min)",gp=gpar(fontsize=15,font=2),rot=90),
# bottom= textGrob("Temperature (\u00B0C)",gp=gpar(fontsize=15,font=2)))

###############GAM GRAPHS#################
smr.temp<-visreg(SMR_gam.1,by="Pond",xvar="SMR.Temp.true", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("Temperature (\u00B0C)")+ylab("")+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+scale_fill_manual(values=c("#FFC0CB80", "#FFD70080","#ADD8E680"))+theme(legend.position="none",strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
smr.pond<-visreg(SMR_gam.1,xvar="Pond", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("")+ylab(NULL)+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+theme(strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+ylim(0,4)
smr.mass<-visreg(SMR_gam.1,xvar="Mass", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("")+ylab("")+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+theme(strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
grid.arrange(
  smr.pond,smr.mass,smr.temp,
  ncol=3
)

fill_alpha(c("pink","gold","lightblue"),alpha=.5)
mmr.temp<-visreg(mmr.mod_gam.1,by="Pond",xvar="MMR.Temp.true", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("")+ylab("")+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+scale_fill_manual(values=c("#FFC0CB80", "#FFD70080","#ADD8E680"))+theme(legend.position="none",strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
mmr.pond<-visreg(mmr.mod_gam.1,xvar="Pond", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("")+ylab(NULL)+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+theme(strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+ylim(0,6)
mmr.mass<-visreg(mmr.mod_gam.1,xvar="Mass", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("")+ylab("")+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+theme(strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
grid.arrange(
  mmr.pond,mmr.mass,mmr.temp,
  ncol=3
)


as.temp<-visreg(as.mod_gam.1,by="Pond",xvar="AS.Temp.true", alpha=0.05, overlay=TRUE,gg=TRUE)+theme_classic()+xlab("Temperature (\u00B0C)")+ylab("")+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+scale_fill_manual(values=c("#FFC0CB80", "#FFD70080","#ADD8E680"))+theme(legend.position="none",strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
as.pond<-visreg(as.mod_gam.1,xvar="Pond", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("Rearing Pond")+ylab(NULL)+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+theme(strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+ylim(0,4)
as.mass<-visreg(as.mod_gam.1,xvar="Mass", alpha=0.05, overlay=TRUE,gg=TRUE)+theme_classic()+xlab("Mass (g)")+ylab("")+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+theme(strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
grid.arrange(
  as.pond,as.mass,as.temp,
  ncol=3
)

expression(paste('Metabolic Rate '(mg*O[2] * kg^-1 * min^-1)))
grid.arrange(
  smr.pond,smr.graph,
  mmr.pond,mmr.graph,
  as.pond,as.graph,
  ncol=2,
  left = textGrob(expression(paste(bold('Metabolic Rate '(mg*O[2] * min^-1)))),rot=90,gp=gpar(fontsize=15))
)






#####################################
#Acclimated vs Acute Respirometry
#####################################
acclim<-filter(final_data, Trial=="Acclimated")

resp.comp<-filter(final_data,Temp.factor!="27"&Temp.factor!="29")
resp.comp$Trial<-factor(resp.comp$Trial,levels=c("Acute","Acclimated"))
resp.comp$Temp.factor<-factor(resp.comp$Temp.factor,levels=c("21","23","25"))
DB.comp<-filter(resp.comp,Pond=="DB")
BC.comp<-filter(resp.comp,Pond=="BC")
ME.comp<-filter(resp.comp,Pond=="ME")

grid.arrange(
  ggplot(data=DB.comp, aes(x = as.factor(Temp.factor), y = SMR.ma, color=Trial)) + geom_boxplot() + geom_point(position = position_jitterdodge()) +ggtitle("Dearborn, MI")+labs(x=NULL, y="Standard Metabolic Rate (mgO2/kg/min)")+ theme(legend.position="none",plot.title = element_text(hjust=.5)) + ylim(0,10)+scale_color_manual(values=c("deeppink2","darkred")),
  ggplot(data=BC.comp, aes(x = as.factor(Temp.factor), y = SMR.ma, color=Trial)) + geom_boxplot() + geom_point(position = position_jitterdodge())+ ggtitle("Bay City, MI")+labs(x=NULL, y=NULL) + theme(legend.position="none",plot.title = element_text(hjust=.5)) + ylim(0,10)+scale_color_manual(values=c("orange","darkorange3")),
  ggplot(data=ME.comp, aes(x = as.factor(Temp.factor), y = SMR.ma, color=Trial)) + geom_boxplot() + geom_point(position = position_jitterdodge()) +ggtitle("Manistique, MI")+ labs(x=NULL, y=NULL)+ theme(legend.position="none",plot.title = element_text(hjust=.5)) + ylim(0,10)+scale_color_manual(values=c("deepskyblue","blue4")),
  ggplot(data=DB.comp, aes(x = as.factor(Temp.factor), y = MMR.Best.ma, color=Trial)) + geom_boxplot() + geom_point(position = position_jitterdodge()) + labs(x=NULL, y="Maximum Metabolic Rate (mgO2/kg/min)")+ theme(legend.position="none") + ylim(0,18)+scale_color_manual(values=c("deeppink2","darkred")),
  ggplot(data=BC.comp, aes(x = as.factor(Temp.factor), y = MMR.Best.ma, color=Trial)) + geom_boxplot() + geom_point(position = position_jitterdodge()) + labs(x=NULL, y=NULL) + theme(legend.position="none") + ylim(0,18)+scale_color_manual(values=c("orange","darkorange3")),
  ggplot(data=ME.comp, aes(x = as.factor(Temp.factor), y = MMR.Best.ma, color=Trial)) + geom_boxplot() + geom_point(position = position_jitterdodge()) + labs(x=NULL, y=NULL)+ theme(legend.position="none")+ ylim(0,18)+scale_color_manual(values=c("deepskyblue","blue4")),
  ggplot(data=DB.comp, aes(x = as.factor(Temp.factor), y = AS.ma, color=Trial)) + geom_boxplot() + geom_point(position = position_jitterdodge()) + labs(x="Temperature, \u00B0C", y="Aerobic Scope (mgO2/kg/min)") + ylim(0,13)+scale_color_manual(values=c("deeppink2","darkred"))+ theme(legend.position="none"),
  ggplot(data=BC.comp, aes(x = as.factor(Temp.factor), y = AS.ma, color=Trial)) + geom_boxplot() + geom_point(position = position_jitterdodge())  + labs(x="Temperature, \u00B0C", y=NULL) + ylim(0,13)+scale_color_manual(values=c("orange","darkorange3"))+ theme(legend.position="none"),
  ggplot(data=ME.comp, aes(x = as.factor(Temp.factor), y = AS.ma, color=Trial)) + geom_boxplot() + geom_point(position = position_jitterdodge())  + labs(x="Temperature, \u00B0C", y=NULL) + ylim(0,13)+scale_color_manual(values=c("deepskyblue","blue4"))+ theme(legend.position="none"),
  ncol=3
)

#####################################
#SMR
#####################################

#Just Acclimated
smr.acclim<-gls(log(SMR.abs)~log(Mass)*Temp.factor*Pond,data=acclim)

sel<-dredge(smr.acclim,rank="AIC",fixed=c("log(Mass)","Temp.factor"),trace=2)
get.models(sel,subset=delta<3)

smr.acclim.final<-gls(log(SMR.abs)~1+Pond+Temp.factor+log(Mass),data=acclim)
anova(smr.acclim.final)
summary(smr.acclim.final)
tab_model(smr.acclim.final)




emms<-emmip(smr.acclim.final,Pond~Temp.factor|log(Mass),plotit=F,CIs=T,style="factor")
smr.emms<-emmip_ggplot(emms)
smr.emms<-smr.emms+theme_classic()+theme(legend.position="None",strip.text.x = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"))+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue","darkred","darkorange3","blue4"))
smr.emms


#Acute vs acclimated acclimated

#Select SMR model through MuMIn

#Global model with all interactions
smr.lme.sel<-lme(log(SMR.abs)~log(Mass)*Temp.factor*Pond*Trial,random=~1|Ind,data=resp.comp)

#SMR model selection
sel.smr<-dredge(smr.lme.sel,rank="AIC",fixed=c("log(Mass)","Temp.factor","Trial","Pond"),trace=2)

#Check all models within delta 2
get.models(sel.smr,subset=delta<2)[[1]]


smr.lme.final<-lme(log(SMR.abs)~Pond + Trial + Temp.factor:Trial + 1 + log(Mass) +      Temp.factor,random=~1|Ind,data=resp.comp)
summary(smr.lme.final)
anova(smr.lme.final)
AIC(smr.lme.final)
tab_model(smr.lme.final)


smr.emm.s<-emmeans(smr.lme.final,~Pond|Temp.factor+Trial+log(Mass))
pairs(smr.emm.s)

smr.emm.s<-emmeans(smr.lme.final,~Trial|Pond+Temp.factor+log(Mass))
pairs(smr.emm.s)


smr.emm.s<-emmeans(smr.lme.final,~Trial+Temp.factor|Pond+log(Mass))
pairs(smr.emm.s)

emms<-emmip(smr.lme.final,Pond+Trial~Temp.factor|log(Mass),plotit=F,CIs=T,type="response")
smr.emms<-emmip_ggplot(emms)
smr.emms<-smr.emms+theme_classic()+theme(legend.position="None",strip.text.x = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+ ggtitle("SMR")+ labs(x="",y="")+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue","darkred","darkorange3","blue4"))
smr.emms


#####################################
#MMR
#####################################

#Just Acclimated
mmr.acclim<-gls(log(MR.Best.abs)~log(Mass)*Temp.factor*Pond,data=acclim)
anova(mmr.acclim)
summary(mmr.acclim)
tab_model(mmr.acclim)

ggplot(acclim,aes(x=Temp.factor,y=MR.Best.abs,color=Pond))+geom_boxplot()


sel<-dredge(mmr.acclim,rank="AIC",fixed=c("log(Mass)","Temp.factor","Pond"),trace=2)
get.models(sel,subset=delta<3)

mmr.acclim.final<-gls(log(MR.Best.abs)~ log(Mass) + Pond + Temp.factor ,data=acclim)
anova(mmr.acclim.final)
summary(mmr.acclim.final)
tab_model(mmr.acclim.final)


emms<-emmip(mmr.acclim.final,Pond~Temp.factor|log(Mass),plotit=F,CIs=T)
mmr.acclim.emms<-emmip_ggplot(emms)
mmr.acclim.emms<-mmr.acclim.emms+theme_classic()+theme(legend.position="None",strip.text.x = element_blank(),strip.text.y=element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"))+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue","darkred","darkorange3","blue4"))
mmr.acclim.emms


#Acute vs acclimated

#Using MuMin to select model
mmr.lme.sel<-lme(log(MR.Best.abs)~log(Mass)*Temp.factor*Trial*Pond,random=~1|Ind,data=resp.comp)
anova(mmr.lme.sel)
summary(mmr.lme.sel)

sel.mmr<-dredge(mmr.lme.sel,rank="AIC",fixed=c("Trial"),trace=2)
get.models(sel.mmr,subset=delta<4)

mmr.lme.final<-lme(log(MR.Best.abs) ~ log(Mass) + Pond + Temp.factor + 1 + Trial ,random=~1|Ind,data=resp.comp)

anova(mmr.lme.final)
summary(mmr.lme.final)
AIC(mmr.lme.final,mmr.lme.final.2,mmr.lme.final.3)
tab_model(mmr.lme.final)

mmr.comp<-plot_model(mmr.lme.final,type="pred",terms=c("Temp.factor","Pond","Mass [10]","Trial"),trace=2)
mmr.comp

mmr.emm.s<-emmeans(mmr.lme.final,~Trial+Pond|Temp.factor+log(Mass))
pairs(mmr.emm.s)


emms<-emmip(mmr.lme.final,Pond+Trial~Temp.factor|log(Mass),plotit = F,CIs=T,type="response")
mmr.emms<-emmip_ggplot(emms)
mmr.emms<-mmr.emms+theme_classic()+theme(strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5),legend.position="none")+ labs(x="",y=expression(paste('Metabolic Rate '(mg*O[2] * min^-1)))) + ggtitle("MMR") +
  scale_color_manual(values=c("deeppink2","orange","deepskyblue","darkred","darkorange3","blue4"))
mmr.emms

#####################################
#AS
#####################################

#Just Acclimated
as.acclim<-gls(log(AS.Best.abs)~log(Mass)*Temp.factor*Pond,data=acclim)

sel<-dredge(as.acclim,rank="AIC",fixed=c("log(Mass)","Temp.factor","Pond"),trace=2)
get.models(sel,subset=delta<3)


as.acclim.final<-gls(log(AS.Best.abs)~1+log(Mass)+Temp.factor+Pond,data=acclim)
anova(as.acclim.final)
summary(as.acclim.final)
AIC(as.acclim.final)
tab_model(as.acclim.final)



emms<-emmip(as.acclim.final,Pond~Temp.factor|log(Mass),plotit=F,CIs=T,type="response")
as.acclim.emms<-emmip_ggplot(emms)
as.acclim.emms<-as.acclim.emms+theme_classic()+theme(legend.position="None",strip.text.x = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"))+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue","darkred","darkorange3","blue4"))
as.acclim.emms


#Acute vs acclimated

#Select model based on MuMin
as.lme.sel<-lme(log(AS.Best.abs)~log(Mass)+Pond*Trial*Temp.factor,random=~1|Ind,data=resp.comp)
anova(as.lme.sel)
summary(as.lme.sel)


sel.as<-dredge(as.lme.sel,rank="AIC",fixed=c("log(Mass)","Temp.factor","Trial","Pond"),trace=2)
get.models(sel.as,subset=delta<2)

as.lme.final<-lme(log(AS.Best.abs) ~1 + log(Mass) + Pond + Temp.factor + Trial, random=~1|Ind,data=resp.comp)

anova(as.lme.final)
AIC(as.lme.final)
tab_model(as.lme.final)


as.emm.s<-emmeans(as.lme.final,~Pond+Trial|Temp.factor+log(Mass))
pairs(as.emm.s)


emms<-emmip(as.lme.final,Pond+Trial~Temp.factor|log(Mass),plotit = F,CIs=T,type="response")
as.emms<-emmip_ggplot(emms)
as.emms<-as.emms+theme_classic()+theme(legend.position="None",strip.text.x = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+labs(x="Temperature (\u00B0C)",y="")+ ggtitle("AS")+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue","darkred","darkorange3","blue4"))
as.emms



################Graphs###################

grid.arrange(
  smr.emms,
  mmr.emms,
  as.emms,
  ncol=1
)


##############Polynomial Acute Models########################
#These were investigated but were not used in our results or discussion in the manuscript

library(ggeffects)

smr.lme.acute<-lme(SMR.abs~poly(Mass,2)+poly(SMR.Temp.true,2)*Pond,random=~1|Ind,data=acute_spring)
anova(smr.lme.acute)
summary(smr.lme.acute)
tab_model(smr.lme.acute)
plot_model(smr.lme.acute)
plot_model(smr.lme.acute,type="pred",terms="Mass [all]")
plot_model(smr.lme.acute,type="pred",terms="SMR.Temp.true [all]")
plot_model(smr.lme.acute,type="pred",terms="Pond [all]")
smr<-plot_model(smr.lme.acute,type="pred",terms=c("SMR.Temp.true [all]","Pond","Mass [7]"))
smr
visreg(SMR_gam.1,scale="response",by="Pond",xvar="SMR.Temp.true",overlay=TRUE,line=list(col=c("deeppink2","orange","deepskyblue")))


RGSMR<-ref_grid(smr.lme.acute,at=list(SMR.Temp.true=c(21,23,25,27,29)))


smr.emm<-emmeans(RGSMR,pairwise~Pond|SMR.Temp.true)
smr.emm


test_predictions(smr.lme.acute,c("SMR.Temp.true","Pond"),p_adjust="tukey")
test_predictions(smr.lme.acute,c("Pond"),p_adjust="tukey")

grid.arrange(
plot_model(smr.lme.acute,type="pred",terms=c("Pond","Mass [10]","SMR.Temp.true[21]")),
plot_model(smr.lme.acute,type="pred",terms=c("Pond","Mass [10]","SMR.Temp.true[23]")),
plot_model(smr.lme.acute,type="pred",terms=c("Pond","Mass [10]","SMR.Temp.true[25]")),
plot_model(smr.lme.acute,type="pred",terms=c("Pond","Mass [10]","SMR.Temp.true[27]")),
plot_model(smr.lme.acute,type="pred",terms=c("Pond","Mass [10]","SMR.Temp.true[29]")),
ncol=3)


smr.21<-predict_response(smr.lme.acute,terms=c("Pond [all]","SMR.Temp.true [21]","Mass [10]"))
test_predictions(smr.21)
smr.23<-predict_response(smr.lme.acute,terms=c("Pond [all]","SMR.Temp.true [23]","Mass [10]"))
test_predictions(smr.23)
smr.25<-predict_response(smr.lme.acute,terms=c("Pond [all]","SMR.Temp.true [25]","Mass [10]"))
test_predictions(smr.25)
smr.27<-predict_response(smr.lme.acute,terms=c("Pond [all]","SMR.Temp.true [27]","Mass [10]"))
test_predictions(smr.27)
smr.29<-predict_response(smr.lme.acute,terms=c("Pond [all]","SMR.Temp.true [29]","Mass [10]"))
test_predictions(smr.29)



test_predictions(smr.lme.acute,c("SMR.Temp.true","Pond"),p_adjust="tukey")
test_predictions(smr.lme.acute,c("Pond"),p_adjust="tukey")




mmr.lme.acute<-lme(MR.Best.abs~poly(Mass,2)+poly(MMR.Temp.true,2)*Pond,random=~1|Ind,data=acute_spring)
anova(mmr.lme.acute)
summary(mmr.lme.acute)
tab_model(mmr.lme.acute)

mmr.lme.acute<-lme(MR.Best.abs~poly(Mass,2)+poly(MMR.Temp.true,2)*Pond+SMR.abs,random=~1|Ind,data=acute_spring)
anova(mmr.lme.acute)
summary(mmr.lme.acute)
tab_model(mmr.lme.acute)

AIC(mmr.lme.acute)

mmr<-plot_model(mmr.lme.acute,type="pred",terms=c("MMR.Temp.true [all]","Pond","Mass[7]"))
mmr
mmr.emm<-emmeans(mmr.lme.acute,pairwise~Pond*MMR.Temp.true|Mass)
mmr.emm
contrast(mmr.emm)
joint_tests(mmr.lme.acute)


RGMMR<-ref_grid(mmr.lme.acute,at=list(MMR.Temp.true=c(21,23,25,27,29)))

mmr.emm<-emmeans(RGMMR,pairwise~Pond|MMR.Temp.true)
mmr.emm

mmr.emm<-emmeans(mmr.lme.acute,pairwise~Pond|MMR.Temp.true+Mass)
mmr.emm

grid.arrange(
plot_model(mmr.lme.acute,type="pred",terms=c("Pond","Mass [10]","MMR.Temp.true[21]")),
plot_model(mmr.lme.acute,type="pred",terms=c("Pond","Mass [10]","MMR.Temp.true[23]")),
plot_model(mmr.lme.acute,type="pred",terms=c("Pond","Mass [10]","MMR.Temp.true[25]")),
plot_model(mmr.lme.acute,type="pred",terms=c("Pond","Mass [10]","MMR.Temp.true[27]")),
plot_model(mmr.lme.acute,type="pred",terms=c("Pond","Mass [10]","MMR.Temp.true[29]")),
ncol=3)



as.lme.acute<-lme(AS.Best.abs~poly(Mass,2)+poly(AS.Temp.true,2)*Pond,random=~1|Ind,data=acute_spring)
anova(as.lme.acute)
summary(as.lme.acute)
tab_model(as.lme.acute)
plot_model(as.lme.acute)
as<-plot_model(as.lme.acute,type="pred",terms=c("AS.Temp.true [all]","Pond","Mass [7]"))
as
visreg(as.mod_gam.1,scale="response",by="Pond",xvar="AS.Temp.true",overlay=TRUE, gg=TRUE)+theme_classic()+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))


as.emm<-emmeans(as.lme.acute,pairwise~Pond|Mass+AS.Temp.true)
as.emm


RGAS<-ref_grid(as.lme.acute,at=list(AS.Temp.true=c(21,23,25,27,29)))

as.emm<-emmeans(RGAS,pairwise~Pond|AS.Temp.true)
as.emm

as.emm<-emmeans(as.lme.acute,pairwise~Pond|AS.Temp.true+Mass)
as.emm


grid.arrange(
  plot_model(as.lme.acute,type="pred",terms=c("Pond","Mass [10]","AS.Temp.true[21]")),
  plot_model(as.lme.acute,type="pred",terms=c("Pond","Mass [10]","AS.Temp.true[23]")),
  plot_model(as.lme.acute,type="pred",terms=c("Pond","Mass [10]","AS.Temp.true[25]")),
  plot_model(as.lme.acute,type="pred",terms=c("Pond","Mass [10]","AS.Temp.true[27]")),
  plot_model(as.lme.acute,type="pred",terms=c("Pond","Mass [10]","AS.Temp.true[29]")),
  ncol=3)

plot_model(as.lme.acute,type="pred",terms=c("Pond","Mass [10]"))


grid.arrange(
  smr.graph,
  smr+theme_classic()+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+scale_fill_manual(values=c("deeppink2","orange","deepskyblue"))+scale_x_continuous(breaks=seq(21,29,by=2)),
  mmr.graph,
  mmr+theme_classic()+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+scale_fill_manual(values=c("deeppink2","orange","deepskyblue"))+scale_x_continuous(breaks=seq(21,29,by=2)),
  as.graph,
  as+theme_classic()+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+scale_fill_manual(values=c("deeppink2","orange","deepskyblue"))+scale_x_continuous(breaks=seq(21,29,by=2))
)

p <- ggpredict(smr.lme.acute, c("SMR.Temp.true [all]", "Pond"))
plot(p)
p <- ggpredict(SMR_gam.1, c("SMR.Temp.true [all]", "Pond"))
plot(p)


smr.temp<-visreg(SMR_gam.1,by="Pond",xvar="SMR.Temp.true", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("Temperature (\u00B0C)")+ylab("")+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+scale_fill_manual(values=c("#FFC0CB80", "#FFD70080","#ADD8E680"))+theme(legend.position="none",strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
smr.pond<-visreg(SMR_gam.1,xvar="Pond", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("")+ylab(NULL)+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+theme(strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+ylim(0,4)
smr.mass<-visreg(SMR_gam.1,xvar="Mass", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("")+ylab("")+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+theme(strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
grid.arrange(
  smr.pond,smr.mass,smr.temp,
  ncol=3
)

fill_alpha(c("pink","gold","lightblue"),alpha=.5)
mmr.temp<-visreg(mmr.mod_gam.1,by="Pond",xvar="MMR.Temp.true", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("")+ylab("")+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+scale_fill_manual(values=c("#FFC0CB80", "#FFD70080","#ADD8E680"))+theme(legend.position="none",strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
mmr.pond<-visreg(mmr.mod_gam.1,xvar="Pond", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("")+ylab(NULL)+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+theme(strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+ylim(0,6)
mmr.mass<-visreg(mmr.mod_gam.1,xvar="Mass", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("")+ylab("")+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+theme(strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
grid.arrange(
  mmr.pond,mmr.mass,mmr.temp,
  ncol=3
)


as.temp<-visreg(as.mod_gam.1,by="Pond",xvar="AS.Temp.true", alpha=0.05, overlay=TRUE,gg=TRUE)+theme_classic()+xlab("Temperature (\u00B0C)")+ylab("")+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+scale_fill_manual(values=c("#FFC0CB80", "#FFD70080","#ADD8E680"))+theme(legend.position="none",strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
as.pond<-visreg(as.mod_gam.1,xvar="Pond", overlay=TRUE,gg=TRUE)+theme_classic()+xlab("Rearing Pond")+ylab(NULL)+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+theme(strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+ylim(0,4)
as.mass<-visreg(as.mod_gam.1,xvar="Mass", alpha=0.05, overlay=TRUE,gg=TRUE)+theme_classic()+xlab("Mass (g)")+ylab("")+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+theme(strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
grid.arrange(
  as.pond,as.mass,as.temp,
  ncol=3
)

expression(paste('Metabolic Rate '(mg*O[2] * kg^-1 * min^-1)))
grid.arrange(
  smr.pond,smr.graph,
  mmr.pond,mmr.graph,
  as.pond,as.graph,
  ncol=2,
  left = textGrob(expression(paste(bold('Metabolic Rate '(mg*O[2] * min^-1)))),rot=90,gp=gpar(fontsize=15))
)



#########Random effect considerations##################
SMR_gam.1<-gam(SMR.abs~Pond+s(Mass,k=5)+s(SMR.Temp.true,by=Pond,k=5)+s(Ind,bs="re"), data=acute_spring, method="REML",select=T)
summary(SMR_gam.1)
SMR_gam.re.1<-gam(SMR.abs~Stock+s(Mass,k=5)+s(SMR.Temp.true,by=Stock,k=5)+s(Pond,bs="re")+s(Ind,bs="re"), data=acute_spring, method="REML",select=T)
summary(SMR_gam.re.1)
SMR_gam.re.2<-gam(SMR.abs~Pond+s(Mass,k=5)+s(SMR.Temp.true,by=Pond,k=5)+s(Stock,bs="re")+s(Ind,bs="re"), data=acute_spring, method="REML",select=T)
summary(SMR_gam.re.2)
acute_spring2<- transform(acute_spring, Pond_in_Stock = interaction(Pond, Stock, drop = TRUE))
SMR_gam.re.3<-gam(SMR.abs~Pond*Stock+s(Mass,k=5)+s(SMR.Temp.true,by=Pond_in_Stock,k=5)+s(Ind,bs="re"), data=acute_spring2, method="REML",select=T)
summary(SMR_gam.re.3)

smr.lme.acute<-lme(SMR.abs~poly(Mass,2)+poly(SMR.Temp.true,2)*Pond,random=~1|Ind,data=acute_spring)
summary(smr.lme.acute)

smr.lme.acute.re<-lme(SMR.abs~poly(Mass,2)+poly(SMR.Temp.true,2)*Pond,random=~1|Stock/Ind,data=acute_spring)
summary(smr.lme.acute.re)
smr.lme.acute.re<-lme(SMR.abs~poly(Mass,2)+poly(SMR.Temp.true,2),random=~1|Stock/Pond/Ind,data=acute_spring)
summary(smr.lme.acute.re)
tab_model(smr.lme.acute.re)

mmr.mod_gam.1<-gam(MR.Best.abs~ s(Mass)+s(Ind,bs="re")+Pond+s(MMR.Temp.true,by=Pond,k=5),
                   data=acute_spring, method="REML",select=T)

mmr.lme.acute<-lme(MR.Best.abs~poly(Mass,2)+poly(MMR.Temp.true,2)*Pond+SMR.abs,random=~1|Ind,data=acute_spring)




as.mod_gam.1<-gam(AS.Best.abs~Pond*Stock+s(Mass,k=5)+s(AS.Temp.true,by=Pond_in_Stock,k=5)+#allow effect of temp factor to have varying intercepts
                    s(Ind,bs="re"),#random effect for individual fish
                  data=acute_spring2, method="REML",select=T)
summary(as.mod_gam.1)
as.lme.acute<-lme(AS.Best.abs~poly(Mass,2)+poly(AS.Temp.true,2)*Pond,random=~1|Ind,data=acute_spring)

