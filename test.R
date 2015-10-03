CI_comp_test	<- function(glm_model=NULL,survival_model=NULL,mediator="mediator",indep_var="indep_var",G=10^4,robust=FALSE,cum_time="max",seed=12345){
  #Required Libraries
  require(timereg)
  require(mvtnorm)
  set.seed(seed)
  #Extract Parameter Estimates
  ##Time Invariant
  gamma		<- survival_model$gamma
  
  ##Time Variant
  cum			<- survival_model$cum
  cum_covar 	<- survival_model$covariance #Always robust
  
  ##GLM Coeffiecients
  alpha 		<- glm_model$coefficients
  alpha_var	<- summary(glm_model)$coefficients[,"Std. Error"]^2
  
  ##Robust Alternative
  ###Not Robust
  if(!robust){ 
    #Time Invariant
    gamma_var	<- survival_model$var.gamma
    #Time Variant
    cum_var		<- survival_model$var.cum
  }
  ###Yes Robust
  if(robust){ 
    #Time Invariant
    gamma_var	<- survival_model$robvar.gamma
    #Time Variant
    cum_var		<- survival_model$robvar.cum
  }
  
  #Determine Variable Status				
  ##Extract Variable Names
  vars	<- c(names(as.data.frame(cum_var)),names(as.data.frame(gamma_var)))
  vars	<- vars[2:length(vars)]
  
  ##Create Variable Frame
  frame				<- as.data.frame(matrix(as.numeric(0),nrow=length(vars),ncol=5))
  names(frame)		<- c("variable","constant",
                     "indep_var","mediator","control")
  frame[,1]			<- vars						
  frame[,5] 			<- 1
  frame$lambda		<- as.numeric(NA)
  frame$lambda_var	<- as.numeric(NA)
  frame$alpha			<- as.numeric(NA)
  frame$alpha_var		<- as.numeric(NA)
  ##For Loop
  for(i in 1:length(frame[,1])){
    #Time Invariant?
    ##Yes
    if(substr(frame[i,1],1,5)=="const"){
      frame[i,2] <- 1
      if(substr(frame[i,1],7,(nchar(frame[i,1])-1))==mediator){
        frame[i,4] <- 1
        frame[i,5] <- 0
      }
      if(substr(frame[i,1],7,(nchar(frame[i,1])-1))==indep_var){
        frame[i,3] <- 1
        frame[i,5] <- 0
      }
      ##No
    } else{
      if(substr(frame[i,1],1,nchar(frame[i,1]))==mediator){
        frame[i,4] <- 1
        frame[i,5] <- 0
      }
      if(substr(frame[i,1],1,nchar(frame[i,1]))==indep_var){
        frame[i,3] <- 1
        frame[i,5] <- 0
      }
    }
  }
  #Mediator and Independent Variable Time Treatment Identical?
  if(frame$constant[frame$variable==indep_var]==frame$constant[frame$variable==mediator]){
    #Extract Correct Estimate
    ##Lambda
    ###Time Invariant?
    Time_Invariant <- frame$constant[frame$variable==indep_var]
    ####Yes
    if(max(frame$constant)==1){
      frame$lambda[frame$constant ==1 ]		<- gamma
      frame$lambda_var[frame$constant==1] 	<- diag(gamma_var)
    }
    ####No
    if(min(frame$constant)==0){
      #Time snapshot for time varying effect
      ##Default; Take Last Value
      #if(cum_time=="median"){
       # require(matrixStats)
       # frame$lambda[frame$constant==0]		<- colMedians(cum[,2:ncol(cum)],na.rm=TRUE)
       # for(i in 1:ncol(cum))
        #row	<- which(abs(as.data.frame(cum)$time-cum_time)==min(abs(as.data.frame(cum)$time-cum_time)))
      # frame$lambda_var[frame$constant==0]	<- cum_var[row,2:ncol(cum_var)]
      #  cum_covar			<- as.data.frame(cum_covar[[row]])
      #}else{
      if(cum_time=="max"){
        row	<- nrow(cum)
        ##Nearest User Specified Value
     
        cum[row,2:ncol(cum)]
      }else{
        row	<- which(abs(as.data.frame(cum)$time-cum_time)==min(abs(as.data.frame(cum)$time-cum_time)))
      }
      frame$lambda[frame$constant==0]		<- cum[row,2:ncol(cum)]
      frame$lambda_var[frame$constant==0]	<- cum_var[row,2:ncol(cum_var)]
      cum_covar			<- as.data.frame(cum_covar[[row]])
    }
    ##Covariances for Omega
    ###Time Invariant?
    ####Yes
    if(Time_Invariant==1){
      covar		<- as.data.frame(gamma_var)
      covar$names <- names(covar)
      covar11 <- covar[covar$names==frame$variable[frame$indep_var==1],frame$variable[frame$indep_var==1]]
      covar12 <- covar[covar$names==frame$variable[frame$indep_var==1],frame$variable[frame$mediator==1]]
      covar22 <- covar[covar$names==frame$variable[frame$mediator==1],frame$variable[frame$mediator==1]]
      ####No
    }else{
      names(cum_covar)	<- names(as.data.frame(cum_var))[2:ncol(cum_var)]
      cum_covar$names		<- names(as.data.frame(cum_var))[2:ncol(cum_var)]
      covar11				<- cum_covar[cum_covar$names==frame$variable[frame$indep_var==1],frame$variable[frame$indep_var==1]]
      covar12				<- cum_covar[cum_covar$names==frame$variable[frame$indep_var==1],frame$variable[frame$mediator==1]]
      covar22				<- cum_covar[cum_covar$names==frame$variable[frame$mediator==1],frame$variable[frame$mediator==1]]
    }
    Omega <- matrix(c(	covar11,covar12,
                       covar12,covar22),byrow=TRUE,ncol=2,nrow=2)
    ##Alpha
    frame$alpha[!frame$mediator==1]		<- alpha
    frame$alpha_var[!frame$mediator==1]	<- alpha_var
    #Simulation
    IE	<- rep(0,G);TE <- rep(0,G)
    Q <- rep(0,G) 
    for(i in 1:G){
      lambda <- rmvnorm(1, mean =c(
        frame$lambda[frame$indep_var==1],	#mean_lambda1/ Direct Effect
        frame$lambda[frame$mediator==1]),	#mean_lambda3/ Indirect Effect
        sigma = Omega)
      alpha_ie <- rnorm(1,mean = frame$alpha[frame$indep_var==1], sd = sqrt(frame$alpha_var[frame$indep_var==1])) 
      IE[i] <- lambda[2] * alpha_ie
      TE[i] <- IE[i] + lambda[1] 
      Q[i] <- IE[i]/TE[i]
    }
    newList <- list("Input" = list(	"glm_model" = 		glm_model,
                                    "survival_model" = 	survival_model, "mediator" = 		mediator,
                                    "indep_var" = 		indep_var,
                                    "G" = 				G,
                                    "robust" = 			robust,
                                    "cum_time" =		cum_time),
                    "Extracted Values" = list("parameters" = frame,
                                              "omega" = Omega),
                    "Results" = list(	"Indirect Effect" = list( 
                      "Mean" = mean(IE),
                      "CI" = quantile(IE,c(0.025, 0.975))), 
                      "Total Effect" = list(
                        "Mean" = mean(TE),
                        "CI" = quantile(TE,c(0.025, 0.975))), 
                      "Indirect / Total effect" = list(
                        "Mean" = mean(Q),
                        "CI" = quantile(Q,c(0.025, 0.975))))) 
    print("Indirect Effect") 
    print(mean(IE)) 
    print(quantile(IE,c(0.025, 0.975))) 
    print("Total Effect") 
    print(mean(TE)) 
    print(quantile(TE,c(0.025, 0.975))) 
    print("Indirect / Total effect") 
    print(mean(Q)) 
    print(quantile(Q,c(0.025, 0.975)))
    return(invisible(newList))
  }else{
    print("Your mediator and independent variable treat time differently. The current version of this function does not handle such a situation. Sorry :/")
  }
}
