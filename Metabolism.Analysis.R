pkgs <- c("here", "readr", "janitor", "mgcv", "gratia", "dplyr", "ggplot2",
          "ggrepel","visreg","tibble","mgcViz","tidymv","tidyverse", "nlme", "emmeans", "kableExtra", "GGally", 
          "qqplotr","grid","gridExtra","car","RColorBrewer","r2glmm","lmerTest","ggpubr","multcomp","jtools","lme4","sjPlot","MuMIn")

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

all_23_acute<-filter(final_data,Temp.factor=="23"&Trial=="Acute")
all_23_acclimated<-filter(final_data,Temp.factor=="23"&Trial=="Acclimated")

all_25_acute<-filter(final_data,Temp.factor=="25"&Trial=="Acute")
all_25_acclimated<-filter(final_data,Temp.factor=="25"&Trial=="Acclimated")


mass.smr.21<-lm(log(SMR.abs)~log(Mass),data=all_21)
summary(mass.smr.21)
mass.smr.23<-lm(log(SMR.abs)~log(Mass),data=all_23)
summary(mass.smr.23)
mass.smr.25<-lm(log(SMR.abs)~log(Mass),data=all_25)
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


for (i in 1:nrow(final_data)){
  if(final_data$Pond[i]=="DB"){
    final_data$SMR.ma[i]<-final_data$SMR.units[i]*(final_data$Mass[i]/10)**(1-0.70681)
  } else if (final_data$Pond[i]=="BC") {
    final_data$SMR.ma[i]<-final_data$SMR.units[i]*(final_data$Mass[i]/10)**(1-0.74080)
  } else if (final_data$Pond[i]=="ME") {
    final_data$SMR.ma[i]<-final_data$SMR.units[i]*(final_data$Mass[i]/10)**(1-0.81857)
  }}

for (i in 1:nrow(all_21)){
  if(all_21$Pond[i]=="DB"){
    all_21$SMR.ma[i]<-all_21$SMR.units[i]*(all_21$Mass[i]/10)**(1-.89)
  } else if (all_21$Pond[i]=="BC") {
    all_21$SMR.ma[i]<-all_21$SMR.units[i]*(all_21$Mass[i]/10)**(1-.89)
  } else if (all_21$Pond[i]=="ME") {
    all_21$SMR.ma[i]<-all_21$SMR.units[i]*(all_21$Mass[i]/10)**(1-.89)
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
ggplot(data=all_21.final, aes(log(Mass),y=log(MMR.Best.ma),color=Pond))+geom_point(size=2)+labs(x="Mass",y="SMR")+geom_smooth(method="lm")+theme_classic()+theme(legend.position = "None",axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))

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
ggplot(data=all_21.final, aes(log(Mass),y=log(AS.ma),color=Trial))+geom_point(size=2)+labs(x="Mass",y="SMR")+geom_smooth(method="lm")+theme_classic()+theme(legend.position = "None",axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))

final_data$AS.fact<-(final_data$MR.Best.abs)/final_data$SMR.abs

#Relationship of SMR with MMR

smrvmmr.lm<-lm(log(MMR.abs)~log(SMR.abs)*Temp.factor*log(Mass),data=final_data)
summary(smrvmmr.lm)

smrvmmr.pond<-ggplot(final_data, aes(x = SMR.abs, y = MR.Best.abs, color=Pond)) + 
  theme_classic()+ geom_point() + scale_y_continuous(expand = expansion(mult = c(0, 0)),limits=c(0,20)) + 
  labs(y="MMR (mgO2/min)",x="SMR (mgO2/min)",title="SMR Relationship to MMR") +
  theme(legend.box.background = element_rect(colour = "black"),axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),plot.title = element_text(size=25,face="bold",hjust = 0.5))+
  guides(color=guide_legend(title="Population"))+ scale_color_manual(values=c("deepskyblue","orange","deeppink2"),labels=c('Dearborn','Bay City','Menominee'))+
  geom_smooth(method='lm')
smrvmmr.pond


smrvmmr.temp<-ggplot(final_data, aes(x = SMR.ma, y = MMR.Best.ma, color=Temp.factor)) + 
  theme_classic()+ geom_point() + scale_y_continuous(expand = expansion(mult = c(0, 0)),limits=c(0,20)) + 
  labs(y="",x=NULL) +
  theme(legend.box.background = element_rect(colour = "black"),axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),plot.title = element_text(size=25,face="bold",hjust = 0.5))+
  guides(color=guide_legend(title="Temperature"),shape=guide_legend(title="Population")) + scale_color_manual(values=c("blue","lightblue","gold","darkgoldenrod1","darkred"))+
  scale_shape_manual(values=c(16,17,15),labels=c('Dearborn','Bay City','Menominee'))+geom_smooth(method="gam")

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
  Pond=factor(final_data$Pond,levels=c("DB","BC","ME"),ordered=F)
)

acute_spring<-filter(final_data,Trial=="Acute")
acute_spring_all<-filter(all_data,Trial=="Acute")
acclim_spring<-filter(final_data,Trial=="Acclimated")

#########################################
###Standard Metabolic Rate
#########################################

##Investigate non-linear relationships


#SMR_gam.1<-gam(SMR.ma~Pond+s(Mass,by=Pond)+s(SMR.Temp.true,by=Pond,k=5)+s(Ind,bs="re"), data=acute_spring, method="REML",select=T)

SMR_gam.1<-gam(SMR.abs~Pond+s(Mass)+s(SMR.Temp.true,by=Pond,k=5)+s(Ind,bs="re"), data=acute_spring, method="REML",select=T)



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


#Plot of Effect of Populations
smr.pop<-visreg(SMR_gam.1,xvar='Pond',gg=T) + labs(x="",y="f(Pond)") + ggtitle("SMR")+
  theme_classic()+theme(legend.position = c(0.8,.2),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds"))
smr.pop


#Plot Smooths
smr.gamS<-smooth_estimates(SMR_gam.1)
smr.gamS<-filter(smr.gamS,smr.gamS$smooth!="s(Mass)")
smr.gamS<-filter(smr.gamS,smr.gamS$type=="TPRS")
smr.gamS$smooth<-factor(smr.gamS$smooth,levels=c("s(SMR.Temp.true):PondDB","s(SMR.Temp.true):PondBC","s(SMR.Temp.true):PondME"))
smr.graph<-ggplot(data=smr.gamS, aes(x=SMR.Temp.true,y=est,color=Pond)) +
  geom_line(size=2) + 
  labs(x="",y="f(Temperature)") + ggtitle("SMR")+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue"),labels=c('Dearborn:LP','Bay City:LP','Menominee:UP')) + 
  theme_classic()+theme(legend.position="",legend.text=element_text(size=rel(1)),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds"))+scale_x_continuous(breaks=c(21,23,25,27,29))
smr.graph


##Check difference between smooths

#First create matrix of values to be plugged into models
pdat<-expand.grid(SMR.Temp.true=seq(21,29,length=400),Mass=seq(min(acute_spring$Mass),max(acute_spring$Mass),length=400),
                  Pond=c('DB',"BC","ME"))
Ind<-seq(1,400)
Ind<-c(Ind,Ind,Ind)
pdat<-cbind(pdat,Ind)

#Using smooth_difff equation do a pairwise comparions among the populations
comp1<-smooth_diff(SMR_gam.2,pdat,'DB','BC','Pond')
comp2<-smooth_diff(SMR_gam.2,pdat,'DB','ME','Pond')
comp3<-smooth_diff(SMR_gam.2,pdat,'BC','ME','Pond')
comp.2<-cbind(SMR.Temp.true=seq(21,29,length=400),
              rbind(comp1,comp2,comp3))
smr.c1<-ggplot(comp.2, aes(x = SMR.Temp.true, y = diff, group = pair)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() + ggtitle("SMR")+
  facet_wrap(~ pair, ncol = 3,labeller="label_both") +
  coord_cartesian(ylim = c(-1,1)) +
  labs(x = "", y = '')+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))
smr.c1




###################################
###Maximum Metabolic Rate
###################################


mmr.mod_gam.1<-gam(MR.Best.abs~ s(SMR.abs,k=5)+s(Mass)+s(Ind,bs="fs")+Pond+s(MMR.Temp.true,by=Pond,k=5),
                   data=acute_spring, method="REML",select=T)

#mmr.mod_gam.2<-gam(MMR.Best.ma~ s(Temp.factor, SMR.abs,id="1",k=5,bs="sz")+s(Mass,by=Pond)+s(Ind,bs="re")+Pond+s(MMR.Temp.true,by=Pond,k=5),
 #                  data=acute_spring, method="REML",select=T)


MMR_LME<-lme(log(MMR.Best.ma)~Pond*Temp.factor*log(Mass),random=~1|Ind,data=acute_spring)
anova(MMR_LME)

ggplot(acute_spring,aes(x=SMR.Temp.true,y=SMR.ma,color=Pond))+geom_point()+geom_smooth(method="gam")
ggplot(acute_spring,aes(x=Temp.factor,y=log(MMR.Best.ma),color=Pond))+geom_point()+geom_smooth(method="lm")
ggplot(acute_spring,aes(x=Temp.factor,y=log(AS.ma),color=Pond))+geom_point()



mmr.mod_gam.1
gam.check(mmr.mod_gam.1)#Check the k significance for the model
visreg(mmr.mod_gam.1)
plot(mmr.mod_gam.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(mmr.mod_gam.1)
summary (mmr.mod_gam.1)
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

#Single out effect of population for plot
mmr.pop<-visreg(mmr.mod_gam.1,xvar='Pond',gg=T)+ggtitle("MMR")+labs(x="",y="f(Pond)")+
  theme_classic()+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds")) 
mmr.pop

#Plot Temperature smooths from model
s.gamr<-smooth_estimates(mmr.mod_gam.1)
t.gamr<-s.gamr[!is.na(s.gamr$MMR.Temp.true),]
DB.r<-data.frame(filter(t.gamr,Pond=="DB"))
DB.r<-DB.r[DB.r$est==max(DB.r$est),]
DB.r<-DB.r[,"MMR.Temp.true"]

BC.r<-data.frame(filter(t.gamr,Pond=="BC"))
BC.r<-BC.r[BC.r$est==max(BC.r$est),]
BC.r<-BC.r[,"MMR.Temp.true"]

ME.r<-data.frame(filter(t.gamr,Pond=="ME"))
ME.r<-ME.r[ME.r$est==max(ME.r$est),]
ME.r<-ME.r[,"MMR.Temp.true"]


mmr.graph<-ggplot(data=t.gamr, aes(x=MMR.Temp.true,y=est,color=Pond)) +
  geom_line(size=2,aes(linetype=Pond)) + scale_linetype_manual(values=c("twodash", "solid","solid"))+guides(linetype ="none" )+
  labs(x="",y="f(Temperature)") + ggtitle("MMR")+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue"),labels=c("Dearborn:LP","Bay City:LP","Menominee:UP","Global")) + 
  theme_classic()+theme(legend.position="",legend.text=element_text(size=rel(1)),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Populations"))+
  scale_x_continuous(breaks=c(21,23,25,27,29))
  #geom_vline(xintercept=ME.r, color="deepskyblue",linetype="dashed")
mmr.graph


#Plot SMR smooths from model
s.gamr<-smooth_estimates(mmr.mod_gam.1)
t.gamr<-s.gamr[!is.na(s.gamr$SMR.abs),]
mmr.st.1<-ggplot(data=t.gamr, aes(x=SMR.abs,y=est,color=Temp.factor)) +
  geom_line(size=2) +labs(x="SMR (mgO2/kg/min/\u00B0C)",y="SMR Effect Trend") +
  ggtitle("SMR Effect on MMR")+
  scale_color_manual(values=c("blue","lightblue","gold","darkgoldenrod1","darkred"),labels=c("21\u00B0C","23\u00B0C","25\u00B0C","27\u00B0C","29\u00B0C","Global")) +
  theme_classic()+theme(legend.text=element_text(size=rel(1.5)),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Temperature"))
mmr.st.1


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
mmr.c.1


#######################################
###Aerobic Scope
#######################################

as.mod_gam.1<-gam(AS.Best.abs~Pond+s(Mass)+s(AS.Temp.true,by=Pond,k=5)+#allow effect of temp factor to have varying intercepts
                     s(Ind,bs="re"),#random effect for individual fish
                   data=acute_spring, method="REML",select=T)

AS_LME<-lme(log(AS.ma)~Pond*Temp.factor*log(Mass),random=~1|Ind,data=acute_spring)
anova(AS_LME)

gam.check(as.mod_gam.1)
summary(as.mod_gam.1)
plot(as.mod_gam.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
anova(as.mod_gam.1)
draw(as.mod_gam.1)
appraise(as.mod_gam.1)
visreg(as.mod_gam.1)


#Graph Population Effect
as.pop<-visreg(as.mod_gam.1,xvar='Pond',gg=T)+labs(x="Pond Populations",y="f(Pond)") + ggtitle("AS")+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue")) + 
  theme_classic()+theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds")) 
as.pop

#Graph smooths
s.gama<-smooth_estimates(as.mod_gam.1)
s.gama<-filter(s.gama,s.gama$type=="TPRS")
s.gama<-s.gama[!is.na(s.gama$AS.Temp.true),]
s.gama$smooth<-factor(s.gama$smooth,levels=c("s(AS.Temp.true):PondDB","s(AS.Temp.true):PondBC","s(AS.Temp.true):PondME"))

DB.r<-data.frame(filter(s.gama,Pond=="DB"))
DB.r<-DB.r[DB.r$est==max(DB.r$est),]
DB.r<-DB.r[,"AS.Temp.true"]

BC.r<-data.frame(filter(s.gama,Pond=="BC"))
BC.r<-BC.r[BC.r$est==max(BC.r$est),]
BC.r<-BC.r[,"AS.Temp.true"]

ME.r<-data.frame(filter(s.gama,Pond=="ME"))
ME.r<-ME.r[ME.r$est==max(ME.r$est),]
ME.r<-ME.r[,"AS.Temp.true"]

as.graph<-ggplot(data=s.gama, aes(x=AS.Temp.true,y=est,color=smooth)) +
  geom_line(size=2,aes(linetype=Pond)) + scale_linetype_manual(values=c("solid", "twodash","solid"))+guides(linetype ="none" )+
  labs(x="Temperature (\u00B0C)",y="f(Temperature)") + ggtitle("AS")+scale_x_continuous(breaks=c(21,23,25,27,29))+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue"),labels=c("Dearborn:LP","Bay City:LP","Menominee:UP")) + 
  theme_classic()+theme(legend.position="",legend.text=element_text(size=rel(1)),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Populations"))#+geom_vline(xintercept=ME.r, color="deepskyblue",linetype="dashed")+geom_vline(xintercept=BC.r, color="orange",linetype="dashed")
as.graph

#Check differences among smooths
pdat<-expand.grid(AS.Temp.true=seq(21,29,length=400),Mass=seq(min(acute_spring$Mass),max(acute_spring$Mass),length=400),Pond=c('BC','DB','ME'))
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
  smr.c1,
  mmr.c.1,
  as.comp,
  ncol=1
)




#Acute Metabolic Phenotypes
grid.arrange(
  ggplot(acute_spring, aes(x = as.factor(Temp.factor), y = SMR.ma, color=Pond)) + geom_boxplot(outlier.shape = NA,lwd=.5)  + 
    theme_classic() + scale_y_continuous(expand = expansion(mult = c(0, 0)),limits=c(0,20)) + geom_point(position = position_jitterdodge())+
    labs(y=NULL,x=NULL)+scale_color_manual(values=c("deeppink2","orange","deepskyblue"),labels=c('Dearborn, LP','Bay City, LP','Menominee, UP'))+
    ggtitle("Standard Metabolic Rate")+ theme(legend.position=c(.3,.8),axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Populations")),
  ggplot(acute_spring, aes(x = as.factor(Temp.factor), y = MMR.Best.ma, color=Pond)) +
    geom_boxplot(outlier.shape = NA,lwd=.5)  + theme_classic() + 
    scale_y_continuous(expand = expansion(mult = c(0, 0)),limits=c(0,20)) + 
    geom_point(position = position_jitterdodge())+
    labs(y=NULL,x=NULL)+
    scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+ggtitle("Maximum Metabolic Rate")+
    theme(legend.position="none",axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds")),
  ggplot(acute_spring, aes(x = as.factor(Temp.factor), y = AS.ma, color=Pond)) + 
    geom_boxplot(outlier.shape = NA,lwd=.5)  + theme_classic() + 
    scale_y_continuous(expand = expansion(mult = c(0, 0)),limits=c(0,20)) + 
    geom_point(position = position_jitterdodge())+
    labs(y=NULL,x=NULL)+scale_color_manual(values=c("deeppink2","orange","deepskyblue"))+ggtitle("Aerobic Scope")+
    theme(legend.position="none",axis.text=element_text(size=16),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+guides(color=guide_legend(title="Ponds")),
  ncol=3,
  left= textGrob("Metabolic Rate (mgO2/kg/min)",gp=gpar(fontsize=15,font=2),rot=90),
  bottom= textGrob("Temperature (\u00B0C)",gp=gpar(fontsize=15,font=2)))







#####################################
#Acclimated vs Acute Respirometry
#####################################
acclim<-filter(final_data, Trial=="Acclimated")

resp.comp<-filter(final_data,Temp.factor!="27"&Temp.factor!="29")
resp.comp$Trial<-factor(resp.comp$Trial,levels=c("Acute","Acclimated"))
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
anova(smr.lme.final)
AIC(smr.lme.final)
tab_model(smr.lme.final)


smr.emm.s<-emmeans(smr.lme.final,~Pond|Temp.factor+Trial+log(Mass))
pairs(smr.emm.s)

smr.emm.s<-emmeans(smr.lme.final,~Temp.factor|Pond+Trial+log(Mass))
pairs(smr.emm.s)

emms<-emmip(smr.lme.final,Pond+Trial~Temp.factor|log(Mass),plotit=F,CIs=T)
smr.emms<-emmip_ggplot(emms)
smr.emms<-smr.emms+theme_classic()+theme(legend.position="None",strip.text.x = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+ ggtitle("SMR")+ labs(x="",y="")+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue","darkred","darkorange3","blue4"))
smr.emms


#####################################
#MMR
#####################################

#Just Acclimated
mmr.acclim<-gls(log(MR.Best.abs)~log(Mass)+Temp.factor+Pond+Temp.factor:log(SMR.abs),data=acclim)
anova(mmr.acclim)
summary(mmr.acclim)

sel<-dredge(mmr.acclim,rank="AIC",fixed=c("log(Mass)","Temp.factor","Pond"),trace=2)
get.models(sel,subset=delta<3)

mmr.acclim.final<-gls(log(MR.Best.abs)~log(SMR.abs):Temp.factor + 1 + log(Mass) + Pond + Temp.factor ,data=acclim)
anova(mmr.acclim.final)
summary(mmr.acclim.final)
tab_model(mmr.acclim.final)


emms<-emmip(mmr.acclim.final,Pond~Temp.factor|log(SMR.abs)+log(Mass),plotit=F,CIs=T)
mmr.acclim.emms<-emmip_ggplot(emms)
mmr.acclim.emms<-mmr.acclim.emms+theme_classic()+theme(legend.position="None",strip.text.x = element_blank(),strip.text.y=element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"))+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue","darkred","darkorange3","blue4"))
mmr.acclim.emms


#Acute vs acclimated

#Using MuMin to select model
mmr.lme.sel<-lme(log(MR.Best.abs)~log(Mass)*Temp.factor*Trial*Pond*Temp.factor:log(SMR.abs),random=~1|Ind,data=resp.comp)
anova(mmr.lme.sel)
summary(mmr.lme.sel)

sel.mmr<-dredge(mmr.lme.sel,rank="AIC",fixed=c("log(Mass)","Temp.factor","Pond","Trial"),trace=2)
summary(get.models(sel.mmr,subset=delta<2)[[1]])


mmr.lme.final<-lme(log(MR.Best.abs) ~ log(SMR.abs):Temp.factor + 1 + log(Mass) + Pond + Temp.factor + Trial,random=~1|Ind,data=resp.comp)
anova(mmr.lme.final)
summary(mmr.lme.final)
AIC(mmr.lme.final)
tab_model(mmr.lme.final)

mmr.emm.s<-emmeans(mmr.lme.final,~Trial|Pond+log(SMR.abs)+Temp.factor+log(Mass))
pairs(mmr.emm.s)

emms<-emmip(mmr.lme.final,Pond+Trial~Temp.factor|log(SMR.abs)+log(Mass),plotit = F,CIs=T)
mmr.emms<-emmip_ggplot(emms)
mmr.emms<-mmr.emms+theme_classic()+theme(legend.position="None",strip.text.x = element_blank(),strip.text.y = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+ labs(x="") + ggtitle("MMR") +
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
summary(as.acclim)
AIC(as.acclim)
tab_model(as.acclim.final)

emms<-emmip(as.acclim.final,Pond~Temp.factor|log(Mass),plotit=F,CIs=T)
as.acclim.emms<-emmip_ggplot(emms)
as.acclim.emms<-as.acclim.emms+theme_classic()+theme(legend.position="None",strip.text.x = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"))+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue","darkred","darkorange3","blue4"))
as.acclim.emms


#Acute vs acclimated

#Select model based on MuMin
as.lme.sel<-lme(log(AS.Best.abs)~log(Mass)+Pond+Trial*Temp.factor,random=~1|Ind,data=resp.comp)
anova(as.lme.sel)
summary(as.lme.sel)


sel.as<-dredge(as.lme.sel,rank="AIC",fixed=c("log(Mass)","Temp.factor","Trial","Pond"),trace=2)
summary(get.models(sel.as,subset=1)[[1]])

as.lme.final<-lme(log(AS.Best.abs) ~1 + log(Mass) + Pond + Temp.factor + Trial, random=~1|Ind,data=resp.comp)
anova(as.lme.final)
AIC(as.lme.final)
tab_model(as.lme.final)


as.emm.s<-emmeans(as.lme.final,~Trial|Pond+Temp.factor+log(Mass))
pairs(as.emm.s)


emms<-emmip(as.lme.final,Pond+Trial~Temp.factor|log(Mass),plotit = F,CIs=T)
as.emms<-emmip_ggplot(emms)
as.emms<-as.emms+theme_classic()+theme(legend.position="None",strip.text.x = element_blank(),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title = element_text(size=15,face="bold",hjust = 0.5))+labs(x="Temperature",y="")+ ggtitle("AS")+
  scale_color_manual(values=c("deeppink2","orange","deepskyblue","darkred","darkorange3","blue4"))
as.emms



################Graphs###################

grid.arrange(
  smr.emms,
  mmr.emms,
  as.emms,
  ncol=1
)


