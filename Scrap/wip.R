###Test

### Aalen Mediation

####### ORGANIZATION  #######
file_drive= "C:"
#file_drive="E:"

directory="/Dropbox/Current Projects/DeLongis/Alameda Survival Manuscript/Alameda Data and Initial Code"
directory<-paste0(file_drive,directory)
content_folder = paste0(directory,"/content")
data_folder = paste0(content_folder,"/data")
graphics_folder = paste0(directory,"/graphics")
mplus_folder = paste0(content_folder,"/mplusautomation")
tables_folder = paste0(content_folder,"/tables")
slides_folder = paste0(content_folder,"/slides")

####### LIBRARIES #######


library(beepr)
library(effects)
library(foreign)
library(ggplot2)
library(Hmisc)
library(knitr)
library(lavaan);library(lme4)
library(matrixStats);library(mediation)
library(survival);library(stargazer)
library(timereg)

####### FUNCTIONS #######
source("https://raw.githubusercontent.com/smasongarrison/Functions/master/Mason_Functions.R")

####### OPTIONS #######
options(width=60)
opts_chunk$set(size="small")
#knit_hooks$set(source=listing, output=listing)
#options(warn=-1) #Warnings Off
set.seed(12345)
####### DATA ######
val <- read.csv(paste0(content_folder,"/Data/Alameda_Husband_Wife.csv"))[,2:19]
val <- subset(val, is.na(val$CP_ID)==FALSE)
val$H_Start<-as.numeric(val$H_Start)
val$W_Start<-as.numeric(val$W_Start)
val$W_End<-as.numeric(val$W_End)
val$H_End<-as.numeric(val$H_End)
val$H_Event<-as.numeric(val$H_Event)
val$W_Event<-as.numeric(val$W_Event)

summary(val)
###########################

##test
ols_m <- glm(W_ACT_74 ~ W_ACT_65+H_ACT_65+HH_INCOME_65+W_Start, data=val)

aalen_tv <- aalen(Surv( 
                       time=W_End, event=W_Event) ~ W_ACT_74 + W_ACT_65 +H_ACT_65 + HH_INCOME_65 + W_Start, data = val, robust=T,covariance=1)

glm_model<-ols_m
survival_model<-aalen_tv
#med_var<-"W_ACT_74"
indep_var<-"W_ACT_74"
#indep_var<-"H_ACT_65"
mediator<-"H_ACT_65"
G=10^4
robust=FALSE
cum_time=median(val$W_End[val$W_Event==1],na.rm=TRUE)

a<-CI_comp_test(glm_model=glm_model,survival_model=survival_model,mediator=mediator,indep_var=indep_var,G=10^4,cum_time=cum_time,robust=FALSE)



#########################
CI_comp_test(glm_model=glm_model,survival_model=survival_model,mediator=mediator,indep_var=indep_var,G=10^4,cum_time=cum_time)





Omega <- matrix(c(0.000197^2,-1.12*10^(-9),
                  -1.12*10^(-9),0.000054^2),byrow=TRUE, nrow=2) 


IE <- rep(0,G);TE <- rep(0,G);Q <- rep(0,G) 
for(i in 1:G) {
  lambda <- rmvnorm(1, mean = c(0.000561   #lambda of ses in aalen L1 (indepvar)
                                ,0.000234 #lambda of physical health L3 mediation)
  ), sigma = Omega)
  alpha <- rnorm(1,mean=0.67 # coffecient of ses predicting physical health
                 , sd=0.066) 
  IE[i] <- lambda[2] * alpha
  TE[i] <- IE[i] + lambda[1] 
  Q[i] <- IE[i]/TE[i]
}
print("Indirect Effect") 
print(mean(IE)) 
print(quantile(IE,c(0.025, 0.975))) 
print("Total Effect") 
print(mean(TE)) 
print(quantile(TE,c(0.025, 0.975))) 
print("Indirect / Total effect") 
print(mean(Q)) 
print(quantile(Q,c(0.025, 0.975)))

