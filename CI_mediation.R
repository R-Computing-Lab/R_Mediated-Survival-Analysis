CI_comp<- function(glm_model=NULL,survival_model=NULL,mediator="mediator",indep_var="indep_var",G=10^4,robust=FALSE,cum_time="last",seed=12345, snapshot=TRUE){
###### Set -Up ######  
# Required Libraries
  require(timereg)
  require(mvtnorm)
# Set Seed
  set.seed(seed)
  
###### Extract Parameter Estimates ######
# Time Invariant
  gamma       <- survival_model$gamma
  gamma_var	  <- survival_model$var.gamma
  
# Time Variant
  cum <- survival_model$cum
  cum_var		  <- survival_model$var.cum
  cum_covar 	<- survival_model$covariance #Always robust
  
# GLM Coeffiecients
  alpha 		  <- glm_model$coefficients
  alpha_var	  <- summary(glm_model)$coefficients[,"Std. Error"]^2
  
# Robust Variance?
  if(robust){ 
    gamma_var	<- survival_model$robvar.gamma  # Time Invariant
    cum_var		<- survival_model$robvar.cum} # Time Variant
  
###### Variable Set-Up ######  		
# Extract Variable Names
  vars	      <- c(names(as.data.frame(cum_var)),names(as.data.frame(gamma_var)))
  vars	      <- vars[2:length(vars)]
  
# Create Variable Frame
  frame				     <- as.data.frame(matrix(as.numeric(0),nrow=length(vars),ncol=5))
  names(frame)     <- c("variable","constant",
                     "indep_var","mediator","control")
  frame[,1]			   <- vars						
  frame[,5] 			 <- 1
  frame$lambda		 <- as.numeric(NA)
  frame$lambda_var <- as.numeric(NA)
  frame$alpha			 <- as.numeric(NA)
  frame$alpha_var	 <- as.numeric(NA)
  
#### For Loop ####
  for(i in 1:length(frame[,1])){
  # Time Invariant
    if(substr(frame[i,1],1,5)=="const"){
      frame[i,2] <- 1
  ## Is Mediator?
      if(substr(frame[i,1],7,(nchar(frame[i,1])-1))==mediator){
        frame[i,4] <- 1
        frame[i,5] <- 0
      }
  ## Is Independent Variable?
      if(substr(frame[i,1],7,(nchar(frame[i,1])-1))==indep_var){
        frame[i,3] <- 1
        frame[i,5] <- 0
      }
  # Time Variant
    } else{
  ## Is Mediator?
      if(substr(frame[i,1],1,nchar(frame[i,1]))==mediator){
        frame[i,4] <- 1
        frame[i,5] <- 0
      }
  ## Is Independent Variable?
      if(substr(frame[i,1],1,nchar(frame[i,1]))==indep_var){
        frame[i,3] <- 1
        frame[i,5] <- 0
      }
    }
  }
###### Extract Estimates ######
# Snap shot?
## YES
  if(snapshot){
  # Mediator and Independent Variable Time Treatment Identical?
  ## Yes
    if(frame$constant[frame$variable==indep_var]==frame$constant[frame$variable==mediator]){
    # Extract Correct Estimate
    ## Lambda
    ### Time Invariant?
    Time_Invariant <- frame$constant[frame$variable==indep_var]
    #### Yes
    if(max(frame$constant)==1){
      frame$lambda[frame$constant==1 ]		  <- gamma
      frame$lambda_var[frame$constant==1] 	<- diag(gamma_var)
    }
    #### No
    if(min(frame$constant)==0){
      # Last Value
      if(cum_time=="last"){
        row	<- nrow(cum)
      }else{
      # Nearest User Specified Value
        row	<- which(abs(as.data.frame(cum)$time-cum_time)==min(abs(as.data.frame(cum)$time-cum_time)))
      }
      # Extract Lambda for Given Row
      frame$lambda[frame$constant==0]		<- cum[row,2:ncol(cum)]
      frame$lambda_var[frame$constant==0]	<- cum_var[row,2:ncol(cum_var)]
      cum_covar			<- as.data.frame(cum_covar[[row]])
    }
    ## Covariances for Omega
    ### Time Invariant?
    #### Yes
    if(Time_Invariant==1){
      covar		<- as.data.frame(gamma_var)
      covar$names <- names(covar)
      covar11 <- covar[covar$names==frame$variable[frame$indep_var==1],frame$variable[frame$indep_var==1]]
      covar12 <- covar[covar$names==frame$variable[frame$indep_var==1],frame$variable[frame$mediator==1]]
      covar22 <- covar[covar$names==frame$variable[frame$mediator==1],frame$variable[frame$mediator==1]]
      #### No
    }else{
      names(cum_covar)	<- names(as.data.frame(cum_var))[2:ncol(cum_var)]
      cum_covar$names		<- names(as.data.frame(cum_var))[2:ncol(cum_var)]
      covar11				<- cum_covar[cum_covar$names==frame$variable[frame$indep_var==1],frame$variable[frame$indep_var==1]]
      covar12				<- cum_covar[cum_covar$names==frame$variable[frame$indep_var==1],frame$variable[frame$mediator==1]]
      covar22				<- cum_covar[cum_covar$names==frame$variable[frame$mediator==1],frame$variable[frame$mediator==1]]
    }
    Omega <- matrix(c(covar11,covar12,
                       covar12,covar22),byrow=TRUE,ncol=2,nrow=2)
    ## Alpha
    frame$alpha[!frame$mediator==1]		<- alpha
    frame$alpha_var[!frame$mediator==1]	<- alpha_var
    # Simulation
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
                                    "survival_model" = 	survival_model, 
                                    "mediator" = mediator, 
                                    "indep_var" = indep_var,
                                    "G" = G,
                                    "robust" = robust,
                                    "cum_time" = cum_time),
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
    ## No
    print("Your mediator and independent variable treat time differently. The current version of this function does not handle such a situation. Sorry :/")
  }
}else{
## No
    # Mediator and Independent Variable Time Treatment Identical?
    ## Yes
    if(frame$constant[frame$variable==indep_var]==frame$constant[frame$variable==mediator]){
      # Extract Correct Estimate
      ## Lambda
      ### Time Invariant?
      Time_Invariant <- frame$constant[frame$variable==indep_var]
      #### No
      if(min(frame$constant)==0){
        frame_full<-as.data.frame(cum)
      }
      #### Yes
      if(max(frame$constant)==1){
        frame$lambda[frame$constant==1]		  <- gamma
        l<-frame$lambda[frame$constant==1]
        frame$lambda_var[frame$constant==1] 	<- diag(gamma_var)
        names                                 <-names(as.data.frame(gamma_var))
        
          for(i in 1:length(names)){
            frame_full[names[i]] <- l[i]}
      }
      
      ## Covariances for Omega
      ### Time Invariant?
      #### Yes
      if(Time_Invariant==1){
        covar		<- as.data.frame(gamma_var)
        covar$names <- names(covar)
        covar11 <- covar[covar$names==frame$variable[frame$indep_var==1],frame$variable[frame$indep_var==1]]
        covar12 <- covar[covar$names==frame$variable[frame$indep_var==1],frame$variable[frame$mediator==1]]
        covar22 <- covar[covar$names==frame$variable[frame$mediator==1],frame$variable[frame$mediator==1]]
        #### No
      }
      frameX<-frame
      Results<-data.frame(cum_var[,1])
      names(Results)<-"Time"
      Results$IE<-as.numeric(NA)
      Results$IE_2.5<-as.numeric(NA)
      Results$IE_97.5<-as.numeric(NA)
      Results$TE<-as.numeric(NA)
      Results$TE_2.5<-as.numeric(NA)
      Results$TE_97.5<-as.numeric(NA)
      Results$Q<-as.numeric(NA)
      Results$Q_2.5<-as.numeric(NA)
      Results$Q_97.5<-as.numeric(NA)
      cum_covarz<-cum_covar
      for(i in 1:nrow(cum)){
        cum_covar<-as.data.frame(cum_covarz[[i]])
        names(cum_covar)	<- names(as.data.frame(cum_var))[2:ncol(cum_var)]
        cum_covar$names		<- names(as.data.frame(cum_var))[2:ncol(cum_var)]
        covar11				<- cum_covar[cum_covar$names==frame$variable[frame$indep_var==1],frame$variable[frame$indep_var==1]]
        covar12				<- cum_covar[cum_covar$names==frame$variable[frame$indep_var==1],frame$variable[frame$mediator==1]]
        covar22				<- cum_covar[cum_covar$names==frame$variable[frame$mediator==1],frame$variable[frame$mediator==1]]
        
          
      Omega <- matrix(c(covar11,covar12,
                         covar12,covar22),byrow=TRUE,ncol=2,nrow=2)
      
      ## Alpha
      frame$alpha[!frame$mediator==1]		<- alpha
      frame$alpha_var[!frame$mediator==1]	<- alpha_var
      ## Lamda
      frame$lambda[!frame$constant==1]<-cum[i,2:ncol(cum)]
      frame$lambda_var[!frame$constant==1]<-cum_var[i,2:ncol(cum_var)]
      
      # Simulation
      IE	<- rep(0,G);TE <- rep(0,G)
      Q <- rep(0,G) 
      for(j in 1:G){
        lambda <- rmvnorm(1, mean =c(
          frame$lambda[frame$indep_var==1],	#mean_lambda1/ Direct Effect
          frame$lambda[frame$mediator==1]),	#mean_lambda3/ Indirect Effect
          sigma = Omega)
        alpha_ie <- rnorm(1,mean = frame$alpha[frame$indep_var==1], sd = sqrt(frame$alpha_var[frame$indep_var==1])) 
        IE[j] <- lambda[2] * alpha_ie
        TE[j] <- IE[j] + lambda[1] 
        Q[j] <- IE[j]/TE[j]
      }
      #assign(Omega,paste0("Omega_",i))
     # assign(frame,paste0("frame_",i))

      
      Results$IE[i]<-mean(IE,na.rm=TRUE)
      Results$IE_2.5[i]<-quantile(IE,0.025,na.rm=TRUE)
      Results$IE_97.5[i]<-quantile(IE,0.975,na.rm=TRUE)
      Results$TE[i]<-mean(TE)
      Results$TE_2.5[i]<-quantile(TE,0.025,na.rm=TRUE)
      Results$TE_97.5[i]<-quantile(TE,0.975,na.rm=TRUE)
      Results$Q[i]<-mean(Q,na.rm=TRUE)
      Results$Q_2.5[i]<-quantile(Q,0.025,na.rm=TRUE)
      Results$Q_97.5[i]<-quantile(Q,0.975,na.rm=TRUE)
      }
      return(Results)
    }else{
    ## No
      print("Your mediator and independent variable treat time differently. The current version of this function does not handle such a situation. Sorry :/")
    }
}
}
