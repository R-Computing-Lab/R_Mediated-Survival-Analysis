#######Lange Function ORIGINAL
#
#       
#ols_m <- glm(logphys ~ age + COHA + child + factor(SES), data=ses_dat_m)
#summary(ols_m)
#
#library(timereg)
#aalen_m1 <- aalen(Surv(T2,event) ~ age + COHA + child + logphys + factor(SES), data = ses_dat_m, robust=T )
#aalen_m2 <- aalen(Surv(T2,event) ~ const(age) + const(COHA) + const(child) + const(logphys) + const(factor(SES)), data = ses_dat_m, robust=T )
#CI_comp <- function(mean_lambda1 ,mean_lambda3,covar11,covar12, covar22, mean_alpha , var_alpha, G=10^4)
#{
# require(mvtnorm)
#Omega <- matrix(c(covar11,covar12,covar12,covar22),nrow=2) 
#IE <- rep(0,G);TE <- rep(0,G);Q <- rep(0,G) for(i in 1:G) {
#lambda <- rmvnorm(1, mean = c(mean_lambda1,mean_lambda2), sigma = Omega)
#alpha <- rnorm(1,mean=mean_alpha, sd=sqrt(var_alpha)) 
#IE[i] <- lambda[2] * alpha
#TE[i] <- IE[i] + lambda[1] 
#Q[i] <- IE[i]/TE[i]
#}
#print("IE:") 
#print(mean(IE)) 
#print(quantile(IE,c(0.025, 0.975))) 
#print("TE:") 
#print(mean(TE)) 
#print(quantile(TE,c(0.025, 0.975))) 
#print("Q:") 
#print(mean(Q)) 
#print(quantile(Q,c(0.025, 0.975)))
#}
#f_out(mean_lambda1=0.000561 ,mean_lambda2=0.000234, covar11=0.000197^2, covar12=-1.12*10^(-9), covar22=0.000054^2, mean_alpha=0.67, var_alpha=0.066^2)
###########################

#######Lange Function FIXED
#
#

###Actual Function
CI_comp <- function(glm_model=NULL,survival_model=NULL,med_var="med_var",indep_var="indep_var",G=10^4,med_constant=TRUE,indep_constant=TRUE){
  require(timereg)
  vars<-names(as.data.frame(survival_model$var.gamma))
if(med_constant && indep_constant){
  
  mean_lambda1<-as.data.frame(survival_model$gamma)[paste0("const(",indep_var,")"),] #natural direct effect
  mean_lambda3 <-as.data.frame(survival_model$gamma)[paste0("const(",med_var,")"),] #natural indirect effect
  covar11 <-as.data.frame(survival_model$var.gamma)[paste0("const(",indep_var,")"),paste0("const(",indep_var,")")]
  covar12 <-as.data.frame(survival_model$var.gamma)[paste0("const(",indep_var,")"),paste0("const(",med_var,")")]
  covar22 <-as.data.frame(survival_model$var.gamma)[paste0("const(",med_var,")"),paste0("const(",med_var,")")]
 
  
  mean_alpha <-as.data.frame(glm_model$coefficients)[paste0(med_var),]
  var_alpha<-summary(glm_model)$coefficients[paste0(med_var),"Std. Error"]^2
   require(mvtnorm)
Omega <- matrix(c(covar11,covar12,
                  covar12,covar22),nrow=2) 
IE <- rep(0,G);TE <- rep(0,G);Q <- rep(0,G) 
for(i in 1:G) {
lambda <- rmvnorm(1, mean = c(mean_lambda1,mean_lambda3), sigma = Omega)
alpha <- rnorm(1,mean=mean_alpha, sd=sqrt(var_alpha)) 
IE[i] <- lambda[2] * alpha
TE[i] <- IE[i] + lambda[1] 
Q[i] <- IE[i]/TE[i]
}
print("Natural Indirect Effect") 
print(mean(IE)) 
print(quantile(IE,c(0.025, 0.975))) 
print("Total Effect") 
print(mean(TE)) 
print(quantile(TE,c(0.025, 0.975))) 
print("Indirect / Total effect") 
print(mean(Q)) 
print(quantile(Q,c(0.025, 0.975)))
} 
if(!(med_constant | indep_constant)){ print("Sorry, I haven't coded non-constant effects") }
}
