#Fit main effects non-parametric Dirichlet process mixture model for potential compliance and parametric outcome model (marginal structural model)
#For ENGAGE-type simulated SMART study
library(truncnorm) #Used to draw from truncated normal distribution
library(LaplacesDemon) #Used to draw from inverse-chi-squared distribution
library(MASS) #Used to draw from multivariate non-truncated normal distribution
library(condMVNorm) #Short-hand functions for computing the conditional distribution of a multivariate normally distributed vector.

#################
#Responder mixture membership indicator
computeZresponders <- function(n, H, W_current, D1_current, D2_current, eta_current, beta_D_responders1, beta_D_responders2,X_mat_responders, sigma_current) {
  # computeZresponders draws the latent mixture indicator for each subject from the responder Dirichlet process mixture
  #
  # Args:
  # n: sample size
  # H: number of DP mixture components
  # W_current: current DP mixture weights
  # D1_current (D_{11}) and D2_current (D_{12}) correspond to stage-1 potential compliances 
  # sigma_current: list (H elements) of current pre-truncation covariance matrix estimate for each DP mixture component
  
  #Vector for storing probabilities of being in each mixture component
  Z_prob_h <- rep(0,H)
  
  #Vector for storing the latent mixture component for each subject
  Z_current <- rep(NA,n)
  dens <- matrix(NA,nrow=n,ncol=H)

  for (i in 1:n) {
    for (h in 1:H) {
      dens[i,h] <- dtmvnorm2( cbind(D1_current,
                                    D2_current)[i,],eta_current[h,]+as.vector((c(X_mat_responders[i,]%*%beta_D_responders1,
                                                                                 X_mat_responders[i,]%*%beta_D_responders2))), sigma_current[[h]],c(0,0),c(1,1))
    }
  }
  #Calculate each mixture component indicator
  for (i in 1:n) {#iterate over n individuals
    Z_prob_h <- rep(0,H)
    for (h in 1:H) {#iterate over H mixture components
      
      Z_prob_h[h] <- (dens[i,h]) * (W_current[h])
      
      if (is.infinite(Z_prob_h[h])){
        Z_prob_h[h] = 10
      }
      
      
      
    }
    #If rare numerical problems occur, sample uniformly from the H mixture components
    if (any(is.na(Z_prob_h))|all(Z_prob_h==0)) {
      Z_current[i] <- sample(1:H,size=1,replace=T,prob=c(1,rep(1,H-1)))
    } else {#Otherwise,
      #Sample from categorical distribution the latent mixture component.
      Z_current[i] <- ifelse(any(Z_prob_h<0),1,sample(1:H,1,replace=T,Z_prob_h))
      
    }
  }
  return(Z_current)
}


#################
computeZnonresponders <- function(n,H,W_current,D1_current,D2_current,D5_current,eta_current, beta_D_non_responders41, beta_D_non_responders42, beta_D_non_responders44, X_mat_non_responders4, sigma_current) {
  # computeZnonresponders draws the latent mixture indicator for each subject for the non-responders' Dirichlet process mixture.
  #
  # Args:
  # n: sample size
  # H: number of DP mixture components
  # W_current: current DP mixture weights
  ##  D1_current (D_{11}) and D2_current (D_{12}) correspond to stage-1 potential compliances while D5_current corresponds to the stage-2 potential compliance, either D_{12} or D_{22}
  # sigma_current: list (H elements) of current pre-truncation covariance matrix estimate for each DP mixture component
  
  # Vector for storing probabilities of being in each mixture component
  Z_prob_h <- rep(0,H)
  
  #Vector for storing the latent mixture component for each subject
  Z_current <- rep(NA,n)
  dens <- matrix(NA,nrow=n,ncol=H)

  for (i in 1:n) {
    for (h in 1:H) {
      dens[i,h] <- dtmvnorm2( cbind(D1_current,
                                    D2_current,
                                    D5_current)[i,], eta_current[h,]+as.vector(c(X_mat_non_responders4[i,]%*%beta_D_non_responders41,
                                                                                 X_mat_non_responders4[i,]%*%beta_D_non_responders42,
                                                                                 X_mat_non_responders4[i,]%*%beta_D_non_responders44)), sigma_current[[h]],c(0,0,0),c(1,1,1))
    }
  }
  #Calculate each mixture component indicator
  for (i in 1:n) {#iterate over n individuals
    Z_prob_h <- rep(0,H)
    for (h in 1:H) {#iterate over H mixture components
      
      Z_prob_h[h] <- (dens[i,h]) * (W_current[h])
      
      if (is.infinite(Z_prob_h[h])){
        Z_prob_h[h] = 10
      }
      
      
      
    }
    #If rare numerical problems occur, sample uniformly from the H mixture components
    if (any(is.na(Z_prob_h))|all(Z_prob_h==0)) {
      Z_current[i] <- sample(1:H,size=1,replace=T,prob=c(1,rep(1,H-1)))
    } else {#Otherwise,
      #Sample from categorical distribution the latent mixture component.
      
      Z_current[i] <- ifelse(any(Z_prob_h<0),1,sample(1:H,1,replace=T,Z_prob_h))
      
    }
  }
  return(Z_current)
}


computeWprime <- function(Z_current, H, alpha_current ) {
  #computeWprime computes the w primes
  #
  # Args:
  # Z_current: vector of current mixture component indicators
  # H: number of mixture components
  # alpha_current: current concentration parameter draw
  
  #Create vector to store w primes of length H
  W_prime<-rep(NA,H)
  
  # W_prime for largest integer to 1
  W_prime[H] <- 1
  
  #Create vector which counts the number of observations belonging to each component
  Zequal <- rep(0,H)
  
  #Create vector where hth component is the number of observations which belong to mixture component greater than h
  Zgreater <- rep(0,H)
  
  #Calculate Zequal and Zgreater
  for (h in 1:(H-1)) {
    for (i in 1:length(Z_current)) {
      Zequal[h] <- Zequal[h] + (Z_current[i] == h)
      Zgreater[h] <- Zgreater[h] + (Z_current[i]>h)
    }
    
    Zequal[h] = Zequal[h] + 1
    Zgreater[h] = Zgreater[h]+alpha_current
    
    #Draw w primes
    W_prime[h] <- rbeta(1,Zequal[h],Zgreater[h])
  }
  return(W_prime)
  
}

computeW <- function(W_prime,H) {
  #This function computes W's given W primes and the number of mixture components
  #
  #Args
  # W_prime: from w primes from stick-breaking representation of the Dirichlet process
  # H: number of mixture components
  
  #Returns W
  
  #Compute W
  W <- rep(0,H)
  W[1] <- W_prime[1]
  W_prime[H] <- 1
  for (h in 2:(H)) {
    W[h] <- W_prime[h]*prod(sapply(1:(h-1),function(x) 1-W_prime[x]))
    
    
  }
  return(W)
}

computeAlpha <- function(W_prime,Z_current,alpha_current,H) {
  #Draw concentration parameter using Metropolis-Hastings step
  #
  # Args
  # W_prime: current w primes from stick breaking representation
  # Z_current: mixture component indicator for each individual
  # alpha_current: current draw of the concentration parameter
  # H: number of mixture components
  
  
  alpha_prop <- rgamma(1,1,1)
  alpha_prob_vector <- rep(1,H-1)
  alpha_current1=alpha_current
  
  u <- runif(1,0,1)
  
  for (h in 1:(H-1)) {
    Zequal <- sum(Z_current == h)
    Zgreater = sum(Z_current>h)
    denom <- dbeta(W_prime[h],1+Zequal,alpha_current+Zgreater)
    
    #As long as there aren't numerical problems
    if (denom>0.00000001) {
      alpha_prob_vector[h] <- dbeta(W_prime[h],1+Zequal,alpha_prop+Zgreater)/denom
    } else {
      alpha_prob_vector[h] <- dbeta(W_prime[h],1+Zequal,alpha_prop+Zgreater)
      
    }
  }
  
  alpha_prod <- prod(alpha_prob_vector)
  
  #Take MH step
  if (!is.nan(alpha_prod)) {
    if (u < min(1,alpha_prod)) {
      alpha_current1 <- alpha_prop
    }}
  return(alpha_current1)
}

library(tidyverse)
library(dplyr)
computeBetas <- function(niter,dat,H_responders,H_non_responders) {
  #This is the most important function
  #This function implements the Gibbs sampler
  #It draws missing potential compliances from their posterior using data augmentation
  #and draws the beta coefficients in the outcome model from their posterior (marginal structural model)
  #
  # Args
  # niter: number of MCMC iterations
  # dat: dataset
  # H_responders: number of mixture component for responders
  # H_non_responders: number of mixture components for non-responders
  
  
  ##Extract variables from the dataset
  
  dat <- as.data.frame(dat)
  #Stage-1 assignment (1 or -1)
  a1 <- dat$a1
  
  #Stage-2 assignment (1 or -1)
  a2 <- dat$a2
  
  #Stage 1 response indicator (1 or -1)
  s <- dat$s
  
  #a2[s==1] <- 0
  #Potential compliances
  D1 <- dat$D11
  D2 <- dat$D12
  D3 <- dat$D21
  D4 <- dat$D22
  #y <- log(0.5+dat$outcome)
  y <- dat$y
  #y <- log(0.5+dat$outcome)
  
  #y <- dat$outcome
  #Baseline covariates
  #X1 <- dat$male
  X1 <- dat$GH_baseline
  X2 <-dat$higheduce
  X3 <- dat$smoke_baseline
  #Stage-1 assignment (1 or -1)
  #a1 <- dat[,1]
  
  #Stage-2 assignment (1 or -1)
  #a2 <- dat[,2]
  
  #Stage 1 response indicator (1 or -1)
  #s <- dat[,3]
  
  #Potential compliances
  #D1 <- dat[,4]
  #D2 <- dat[,5]
  #D3 <- dat[,6]
  #D4 <- dat[,7]
  
  #Final outcome
  #y <- dat[,8]
  
  #Baseline covariates
  #X1 <- dat[,10]
  #X2 <-dat[,11]
  #X3 <-dat[,12]
  
  
  #concentration parameters
  alpha_results <- rep(NA,niter)
  
  #Pretruncation mean matrix H mixture components by 3 potential compliances
  eta_current_responders <- matrix(0.5,nrow=H_responders,ncol=2)
  eta_current_non_responders3 <- matrix(0.5,nrow=H_non_responders,ncol=3)
  eta_current_non_responders4 <- matrix(0.5,nrow=H_non_responders,ncol=3)
  
  
  #Potential compliances
  D1_current <- D2_current <- D3_current <- D4_current <-rep(0.5,length(y))
  
  
  #Initialize the membership to each of H mixture components
  Z_current_responders <- sample(1:H_responders,size=length(which(s==1)),replace=T,prob=rep(1,H_responders))
  Z_current_non_responders3 <- sample(1:H_non_responders,size=length(which(s==0&a1==1)),replace=T,prob=rep(1,H_non_responders))
  Z_current_non_responders4 <- sample(1:H_non_responders,size=length(which(s==0&a1==-1)),replace=T,prob=rep(1,H_non_responders))
  
  Z_results <- matrix(NA,nrow=niter,ncol=length(y))
  
  
  #Pre-truncation variances for multivariate truncated normal
  sigma_current_non_responders3 <- vector("list",H_non_responders)
  sigma_current_non_responders4 <- vector("list",H_non_responders)
  sigma_current_responders <- vector("list",H_responders)
  
  for (i in 1:H_responders) {
    sigma_current_responders[[i]] <- diag(1,2)
  }
  for (i in 1:H_non_responders) {
    sigma_current_non_responders3[[i]] <- diag(1,3)
    sigma_current_non_responders4[[i]] <- diag(1,3)
    
  }
  
  #Concentration parameter initialization
  alpha_current_responders <- 1
  alpha_current_non_responders3 <-1
  alpha_current_non_responders4 <-1
  
  #Matrices to store potential compliance imputations
  D1_matrix <-D2_matrix <-D3_matrix <- D4_matrix <- matrix(0.5,nrow=niter,ncol=length(y))
  ############
  
  
  #############
  
  X_seq_overall <- matrix(NA,nrow=length(y),ncol=8)
  
  beta_seq_overall_results <- matrix(NA,nrow=niter,ncol=8)
  beta_seq_overall <-rep(1,8)
  
  #Vector to store residual variance
  sigma_seq_1_squared_results <- rep(NA,niter)
  
  #Residual variance
  sigma_seq_overall_squared <- 1
  
  #Weights for under/overrepresentation in both imputation and fitting marginal structural model.
  w_responders=0.5 #Overrepresented by randomizing once
  w_non_responders=0.25 #underrepresented by randomizing twice
  


  
  
  beta_seq_mean_overall <-beta_seq_mean_overall22 <- MASS::mvrnorm(1,c(0,0,0,0,0,0,0,0),diag(1,8))
  

  beta_D_responders1_sigma_squared<-1
  beta_D_responders2_sigma_squared <-1
  
  beta_D_non_responders44_sigma_squared <- 1
  
  beta_D_non_responders32_sigma_squared <- c(0,0,0)
  beta_D_non_responders33_sigma_squared <- c(0,0,0)
  beta_D_non_responders43_sigma_squared <- c(0,0,0)
  beta_D_non_responders42_sigma_squared <- c(0,0,0)
  beta_D_non_responders41_sigma_squared <- c(0,0,0)
  
  beta_D_non_responders41 <- rep(0.5,3)
  beta_D_non_responders42 <-rep(0.5,3)
  beta_D_non_responders44 <-rep(0.5,3)
  
  beta_D_non_responders31_sigma_squared <-1
  
  beta_D_responders1 <-rep(0.5,3)
  beta_D_responders2 <- rep(0.5,3)
  
  beta_D_non_responders31 <-rep(0.5,3)
  beta_D_non_responders32 <-rep(0.5,3)
  beta_D_non_responders33 <- rep(0.5,3)
  
  
  
  
  
  beta_seq_mean_overall <-beta_seq_mean_overall22 <- MASS::mvrnorm(1,c(0,0,0,0,0,0,0,0),diag(1,8))
  
  X_mat_responders <- cbind(X1,X2,X3)[s==1,]
  X_mat_non_responders3 <- cbind(X1,X2,X3)[s==0&a1==1,]
  X_mat_non_responders4 <- cbind(X1,X2,X3)[s==0&a1==-1,]
  #Iterate through niter MCMC/Gibbs sampler draws
  for (j in 1:niter) {#Each draw imputes the missing potential compliances after fitting their joint posterior

    for (i in 1:length(which(s==1))) {
      #The conditional posteriors are a function of the observed potential compliances and the parametric (normal) outcome (y) models
      
      #Treatment sequence 1
      if (a1[which(s==1)][i]==1) {
        
        #Observed compliance
        D1_current[which(s==1)][i] = D1[which(s==1)][i]
        
        #Conditional mean and variance of full conditional for unobserved compliance, D2
        
        cond_mean <- condMVN(c(X_mat_responders[i,]%*%beta_D_responders1,
                               X_mat_responders[i,]%*%beta_D_responders2)+eta_current_responders[Z_current_responders[i],],sigma_current_responders[[Z_current_responders[i]]],dependent.ind = 2,given.ind = c(1),c(D1_current[which(s==1)][i]),check.sigma = F)$condMean
        cond_var <- condMVN(c(X_mat_responders[i,]%*%beta_D_responders1,
                              X_mat_responders[i,]%*%beta_D_responders2)+eta_current_responders[Z_current_responders[i],],sigma_current_responders[[Z_current_responders[i]]],dependent.ind = 2,given.ind = c(1),c(D1_current[which(s==1)][i]),check.sigma = F)$condVar
        
        
        D2_cond_mean <-(a1[which(s==1)][i]*((beta_seq_overall[3])*(y[which(s==1)][i]-beta_seq_overall[1]-beta_seq_overall[2]*D1_current[which(s==1)][i]*a1[which(s==1)][i]-beta_seq_overall[6]*(X1[which(s==1)][i])-beta_seq_overall[7]*(X2[which(s==1)][i])-beta_seq_overall[8]*(X3[which(s==1)][i]))*cond_var[[1]])+w_responders*sigma_seq_overall_squared*cond_mean)/((beta_seq_overall[3])^2*cond_var+w_responders*sigma_seq_overall_squared)
        D2_cond_var <- w_responders*sigma_seq_overall_squared*cond_var/((beta_seq_overall[3])^2*cond_var+w_responders*sigma_seq_overall_squared)  
        
        
        
        
        if (D2_cond_var <=0) {#Check for numerical problems
          D2_current[which(s==1)][i] <- D2_cond_mean
        } else {#If no numerical problems, impute D2
          D2_current[which(s==1)][i] <- rtruncnorm(1,0,1,D2_cond_mean,sqrt(D2_cond_var))
        }
        
        D3_current[which(s==1)][i] <-0
        
        
        D4_current[which(s==1)][i] <-0
        
        # }
        
        
        
      }
      ######
      
      ##########
      #Treatment sequence 4
      if (a1[which(s==1)][i]==-1) {
        
        #Observed compliance
        D2_current[which(s==1)][i] = D2[which(s==1)][i]
        
        #Conditional mean and variance of full conditional for unobserved compliance
        cond_mean <- condMVN(c(X_mat_responders[i,]%*%beta_D_responders1,
                               X_mat_responders[i,]%*%beta_D_responders2)+eta_current_responders[Z_current_responders[i],],sigma_current_responders[[Z_current_responders[i]]],dependent.ind = 1,given.ind = c(2),c(D2_current[which(s==1)][i]),check.sigma = F)$condMean
        cond_var <- condMVN(c(X_mat_responders[i,]%*%beta_D_responders1,
                              X_mat_responders[i,]%*%beta_D_responders2)+eta_current_responders[Z_current_responders[i],],sigma_current_responders[[Z_current_responders[i]]],dependent.ind = 1,given.ind = c(2),c(D2_current[which(s==1)][i]),check.sigma = F)$condVar
        
        #Conditioning on potential outcome as D1 is in outcome model for this treatment sequence
        D1_cond_mean <-(a1[which(s==1)][i]*((beta_seq_overall[2])*(y[which(s==1)][i]-beta_seq_overall[1]-beta_seq_overall[3]*D2_current[which(s==1)][i]*a1[which(s==1)][i]-beta_seq_overall[6]*(X1[which(s==1)][i])-beta_seq_overall[7]*(X2[which(s==1)][i])-beta_seq_overall[8]*(X3[which(s==1)][i]))*cond_var[[1]])+w_responders*sigma_seq_overall_squared*cond_mean)/((beta_seq_overall[2])^2*cond_var+w_responders*sigma_seq_overall_squared)
        D1_cond_var <- w_responders*sigma_seq_overall_squared*cond_var/((beta_seq_overall[2])^2*cond_var+w_responders*sigma_seq_overall_squared)  
        
        
        if (D1_cond_var <=0) {
          D1_current[which(s==1)][i] <- D1_cond_mean
        } else {
          D1_current[which(s==1)][i] <- rtruncnorm(1,0,1,D1_cond_mean,sqrt(D1_cond_var))
        }
        
        
        D3_current[which(s==1)][i] <-0
        
        D4_current[which(s==1)][i] <-0
        
        
        
      }
      
      
    }
    
    
    
    
    
    for (i in 1:length(which(s==0&a1==1))){
      
      ##########
      #Treatment sequence 2
      if (a1[which(s==0&a1==1)][i]==1&a2[which(s==0&a1==1)][i]==1) {
        
        #Observed compliance
        D1_current[which(s==0&a1==1)][i] = D1[which(s==0&a1==1)][i]
        D3_current[which(s==0&a1==1)][i] = D3[which(s==0&a1==1)][i]
        
        #Conditional mean and variance of full conditional for unobserved compliance, D2
        cond_mean <- condMVN(c((X_mat_non_responders3[i,])%*%beta_D_non_responders31,
                               (X_mat_non_responders3[i,])%*%beta_D_non_responders32,
                               (X_mat_non_responders3[i,])%*%beta_D_non_responders33)+eta_current_non_responders3[Z_current_non_responders3[i],],sigma_current_non_responders3[[Z_current_non_responders3[i]]],dependent.ind = 2,given.ind = c(1,3),c(D1_current[which(s==0&a1==1)][i],D3_current[which(s==0&a1==1)][i]),check.sigma = F)$condMean
        cond_var <- condMVN(c((X_mat_non_responders3[i,])%*%beta_D_non_responders31,
                              (X_mat_non_responders3[i,])%*%beta_D_non_responders32,
                              (X_mat_non_responders3[i,])%*%beta_D_non_responders33)+eta_current_non_responders3[Z_current_non_responders3[i],],sigma_current_non_responders3[[Z_current_non_responders3[i]]],dependent.ind = 2,given.ind = c(1,3),c(D1_current[which(s==0&a1==1)][i],D3_current[which(s==0&a1==1)][i]),check.sigma = F)$condVar
        
        D2_cond_mean <-(a1[which(s==0&a1==1)][i]*((beta_seq_overall[3])*(y[which(s==0&a1==1)][i]-beta_seq_overall[1]-beta_seq_overall[2]*D1_current[which(s==0&a1==1)][i]*a1[which(s==0&a1==1)][i]-beta_seq_overall[4]*D3_current[which(s==0&a1==1)][i]*a2[which(s==0&a1==1)][i]-beta_seq_overall[6]*(X1[which(s==0&a1==1)][i])-beta_seq_overall[7]*(X2[which(s==0&a1==1)][i])-beta_seq_overall[8]*(X3[which(s==0&a1==1)][i]))*cond_var[[1]])+w_non_responders*sigma_seq_overall_squared*cond_mean)/((beta_seq_overall[3])^2*cond_var+w_non_responders*sigma_seq_overall_squared)
        D2_cond_var <- w_non_responders*sigma_seq_overall_squared*cond_var/((beta_seq_overall[3])^2*cond_var+w_non_responders*sigma_seq_overall_squared)  
        
        
        if (D2_cond_var <=0) {
          D2_current[which(s==0&a1==1)][i] <- D2_cond_mean
        } else {
          D2_current[which(s==0&a1==1)][i] <- rtruncnorm(1,0,1,D2_cond_mean,sqrt(D2_cond_var))
        }
        
        
        #if (D4_cond_var <=0) {
        D4_current[which(s==0&a1==1)][i] <- 0
        #} else {
        D4_current[which(s==0&a1==1)][i] <-0
        #}
      }
      ##########
      #Treatment sequence 3
      if (a1[which(s==0&a1==1)][i]==1&a2[which(s==0&a1==1)][i]==-1) {
        #D3 is modeled in the outcome model, so we condition on the outcome, y.
        
        
        #Observed compliance
        D1_current[which(s==0&a1==1)][i] = D1[which(s==0&a1==1)][i]
        
        #Conditional mean and variance of full conditional for unobserved compliance
        
        cond_mean <- condMVN(c((X_mat_non_responders3[i,])%*%beta_D_non_responders31,
                               (X_mat_non_responders3[i,])%*%beta_D_non_responders32,
                               (X_mat_non_responders3[i,])%*%beta_D_non_responders33)+eta_current_non_responders3[Z_current_non_responders3[i],],sigma_current_non_responders3[[Z_current_non_responders3[i]]],dependent.ind = 3,given.ind = c(1,2),c(D1_current[which(s==0&a1==1)][i],D2_current[which(s==0&a1==1)][i]),check.sigma = F)$condMean
        cond_var <-  condMVN(c((X_mat_non_responders3[i,])%*%beta_D_non_responders31,
                               (X_mat_non_responders3[i,])%*%beta_D_non_responders32,
                               (X_mat_non_responders3[i,])%*%beta_D_non_responders33)+eta_current_non_responders3[Z_current_non_responders3[i],],sigma_current_non_responders3[[Z_current_non_responders3[i]]],dependent.ind = 3,given.ind = c(1,2),c(D1_current[which(s==0&a1==1)][i],D2_current[which(s==0&a1==1)][i]),check.sigma = F)$condVar
        
        D3_cond_mean <-(a2[which(s==0&a1==1)][i]*((beta_seq_overall[4])*(y[which(s==0&a1==1)][i]-beta_seq_overall[1]-beta_seq_overall[2]*D1_current[which(s==0&a1==1)][i]*a1[which(s==0&a1==1)][i]-beta_seq_overall[3]*D2_current[which(s==0&a1==1)][i]*a1[which(s==0&a1==1)][i]-beta_seq_overall[6]*(X1[which(s==0&a1==1)][i])-beta_seq_overall[7]*(X2[which(s==0&a1==1)][i])-beta_seq_overall[8]*(X3[which(s==0&a1==1)][i]))*cond_var[[1]])+w_non_responders*sigma_seq_overall_squared*cond_mean)/((beta_seq_overall[4])^2*cond_var+w_non_responders*sigma_seq_overall_squared)
        D3_cond_var <- w_non_responders*sigma_seq_overall_squared*cond_var/((beta_seq_overall[4])^2*cond_var+w_non_responders*sigma_seq_overall_squared)  
        
        
        if (D3_cond_var <=0) {
          D3_current[which(s==0&a1==1)][i] <- D3_cond_mean
        } else {
          D3_current[which(s==0&a1==1)][i] <- rtruncnorm(1,0,1,D3_cond_mean,sqrt(D3_cond_var))
        }
        
        
        
        D4_current[which(s==0&a1==1)][i] <- 0
        
        
        #Conditional mean and variance of full conditional for unobserved compliance, D2
        cond_mean <- condMVN(c((X_mat_non_responders3[i,])%*%beta_D_non_responders31,
                               (X_mat_non_responders3[i,])%*%beta_D_non_responders32,
                               (X_mat_non_responders3[i,])%*%beta_D_non_responders33)+eta_current_non_responders3[Z_current_non_responders3[i],],sigma_current_non_responders3[[Z_current_non_responders3[i]]],dependent.ind = 2,given.ind = c(1,3),c(D1_current[which(s==0&a1==1)][i],D3_current[which(s==0&a1==1)][i]),check.sigma = F)$condMean
        cond_var <-  condMVN(c((X_mat_non_responders3[i,])%*%beta_D_non_responders31,
                               (X_mat_non_responders3[i,])%*%beta_D_non_responders32,
                               (X_mat_non_responders3[i,])%*%beta_D_non_responders33)+eta_current_non_responders3[Z_current_non_responders3[i],],sigma_current_non_responders3[[Z_current_non_responders3[i]]],dependent.ind = 2,given.ind = c(1,3),c(D1_current[which(s==0&a1==1)][i],D3_current[which(s==0&a1==1)][i]),check.sigma = F)$condVar
        
        D2_cond_mean <-(a1[which(s==0&a1==1)][i]*((beta_seq_overall[3])*(y[which(s==0&a1==1)][i]-beta_seq_overall[1]-beta_seq_overall[2]*D1_current[which(s==0&a1==1)][i]*a1[which(s==0&a1==1)][i]-beta_seq_overall[4]*D3_current[which(s==0&a1==1)][i]*a2[which(s==0&a1==1)][i]-beta_seq_overall[6]*(X1[which(s==0&a1==1)][i])-beta_seq_overall[7]*(X2[which(s==0&a1==1)][i])-beta_seq_overall[8]*(X3[which(s==0&a1==1)][i]))*cond_var[[1]])+w_non_responders*sigma_seq_overall_squared*cond_mean)/((beta_seq_overall[3])^2*cond_var+w_non_responders*sigma_seq_overall_squared)
        D2_cond_var <- w_non_responders*sigma_seq_overall_squared*cond_var/((beta_seq_overall[3])^2*cond_var+w_non_responders*sigma_seq_overall_squared)  
        
        
        if (D2_cond_var <=0) {
          D2_current[which(s==0&a1==1)][i] <- D2_cond_mean
        } else {
          D2_current[which(s==0&a1==1)][i] <- rtruncnorm(1,0,1,D2_cond_mean,sqrt(D2_cond_var))
        }
        
        
        
      }
      
    }
    for (i in 1:length(which(s==0&a1==-1))){
      
      ##########
      #Treatment sequence 5
      if (a1[which(s==0&a1==-1)][i]==-1&a2[which(s==0&a1==-1)][i]==1) {
        
        #Observed compliances
        D2_current[which(s==0&a1==-1)][i] = D2[which(s==0&a1==-1)][i]
        D4_current[which(s==0&a1==-1)][i] = D4[which(s==0&a1==-1)][i]
        
        #Conditional mean and variance of full conditional for unobserved compliance, D1
        cond_mean <- condMVN(c((X_mat_non_responders4[i,])%*%beta_D_non_responders41,
                               (X_mat_non_responders4[i,])%*%beta_D_non_responders42,
                               (X_mat_non_responders4[i,])%*%beta_D_non_responders44)+eta_current_non_responders4[Z_current_non_responders4[i],],sigma_current_non_responders4[[Z_current_non_responders4[i]]],dependent.ind = 1,given.ind = c(2,3),c(D2_current[which(s==0&a1==-1)][i],D4_current[which(s==0&a1==-1)][i]),check.sigma = F)$condMean
        cond_var <- condMVN(c((X_mat_non_responders4[i,])%*%beta_D_non_responders41,
                              (X_mat_non_responders4[i,])%*%beta_D_non_responders42,
                              (X_mat_non_responders4[i,])%*%beta_D_non_responders44)+eta_current_non_responders4[Z_current_non_responders4[i],],sigma_current_non_responders4[[Z_current_non_responders4[i]]],dependent.ind = 1,given.ind = c(2,3),c(D2_current[which(s==0&a1==-1)][i],D4_current[which(s==0&a1==-1)][i]),check.sigma = F)$condVar
        
        D1_cond_mean <-(a1[which(s==0&a1==-1)][i]*((beta_seq_overall[2])*(y[which(s==0&a1==-1)][i]-beta_seq_overall[1]-beta_seq_overall[3]*D2_current[which(s==0&a1==-1)][i]*a1[which(s==0&a1==-1)][i]-beta_seq_overall[5]*D4_current[which(s==0&a1==-1)][i]*a2[which(s==0&a1==-1)][i]-beta_seq_overall[6]*(X1[which(s==0&a1==-1)][i])-beta_seq_overall[7]*(X2[which(s==0&a1==-1)][i])-beta_seq_overall[8]*(X3[which(s==0&a1==-1)][i]))*cond_var[[1]])+w_non_responders*sigma_seq_overall_squared*cond_mean)/((beta_seq_overall[2])^2*cond_var+w_non_responders*sigma_seq_overall_squared)
        D1_cond_var <- w_non_responders*sigma_seq_overall_squared*cond_var/((beta_seq_overall[2])^2*cond_var+w_non_responders*sigma_seq_overall_squared)  
        
        
        if (D1_cond_var <=0) {
          D1_current[which(s==0&a1==-1)][i] <- D1_cond_mean
        } else {
          D1_current[which(s==0&a1==-1)][i] <- rtruncnorm(1,0,1,D1_cond_mean,sqrt(D1_cond_var))
        }
        
        
        
        if (D3_cond_var <=0) {
          D3_current[which(s==0&a1==-1)][i] <- 0
        } else {
          D3_current[which(s==0&a1==-1)][i] <- 0
        }
        
      }
      ##########
      #Treatment sequence 6
      if (a1[which(s==0&a1==-1)][i]==-1&a2[which(s==0&a1==-1)][i]==-1) {
        
        #Observed compliances
        D2_current[which(s==0&a1==-1)][i] = D2[which(s==0&a1==-1)][i]
        #Conditional mean and variance of full conditional for unobserved compliance
        cond_mean <- condMVN(c((X_mat_non_responders4[i,])%*%beta_D_non_responders41,
                               (X_mat_non_responders4[i,])%*%beta_D_non_responders42,
                               (X_mat_non_responders4[i,])%*%beta_D_non_responders44)+eta_current_non_responders4[Z_current_non_responders4[i],],sigma_current_non_responders4[[Z_current_non_responders4[i]]],dependent.ind = 3,given.ind = c(1,2),c(D1_current[which(s==0&a1==-1)][i],D2_current[which(s==0&a1==-1)][i]),check.sigma = F)$condMean
        cond_var <- condMVN(c((X_mat_non_responders4[i,])%*%beta_D_non_responders41,
                              (X_mat_non_responders4[i,])%*%beta_D_non_responders42,
                              (X_mat_non_responders4[i,])%*%beta_D_non_responders44)+eta_current_non_responders4[Z_current_non_responders4[i],],sigma_current_non_responders4[[Z_current_non_responders4[i]]],dependent.ind = 3,given.ind = c(1,2),c(D1_current[which(s==0&a1==-1)][i],D2_current[which(s==0&a1==-1)][i]),check.sigma = F)$condVar
        
        
        #Conditioning on potential outcome as D4 is in outcome model for this treatment sequence
        D4_cond_mean <-(a2[which(s==0&a1==-1)][i]*((beta_seq_overall[5])*(y[which(s==0&a1==-1)][i]-beta_seq_overall[1]-beta_seq_overall[2]*D1_current[which(s==0&a1==-1)][i]*a1[which(s==0&a1==-1)][i]-beta_seq_overall[3]*D2_current[which(s==0&a1==-1)][i]*a1[which(s==0&a1==-1)][i]-beta_seq_overall[6]*(X1[which(s==0&a1==-1)][i])-beta_seq_overall[7]*(X2[which(s==0&a1==-1)][i])-beta_seq_overall[8]*(X3[which(s==0&a1==-1)][i]))*cond_var[[1]])+w_non_responders*sigma_seq_overall_squared*cond_mean)/((beta_seq_overall[5])^2*cond_var+w_non_responders*sigma_seq_overall_squared)
        D4_cond_var <- w_non_responders*sigma_seq_overall_squared*cond_var/((beta_seq_overall[5])^2*cond_var+w_non_responders*sigma_seq_overall_squared)  
        
        if (D4_cond_var <=0) {
          D4_current[which(s==0&a1==-1)][i] <- D4_cond_mean
        } else {
          D4_current[which(s==0&a1==-1)][i] <- rtruncnorm(1,0,1,D4_cond_mean,sqrt(D4_cond_var))
        }
        
        
        
        if (cond_var <=0) {
          D3_current[which(s==0&a1==-1)][i] <- 0
        } else {
          D3_current[which(s==0&a1==-1)][i] <- 0
        }
        
        #Conditional mean and variance of full conditional for unobserved compliance, D1
        cond_mean <- condMVN(c((X_mat_non_responders4[i,])%*%beta_D_non_responders41,
                               (X_mat_non_responders4[i,])%*%beta_D_non_responders42,
                               (X_mat_non_responders4[i,])%*%beta_D_non_responders44)+eta_current_non_responders4[Z_current_non_responders4[i],],sigma_current_non_responders4[[Z_current_non_responders4[i]]],dependent.ind = 1,given.ind = c(2,3),c(D2_current[which(s==0&a1==-1)][i],D4_current[which(s==0&a1==-1)][i]),check.sigma = F)$condMean
        cond_var <- condMVN(c((X_mat_non_responders4[i,])%*%beta_D_non_responders41,
                              (X_mat_non_responders4[i,])%*%beta_D_non_responders42,
                              (X_mat_non_responders4[i,])%*%beta_D_non_responders44)+eta_current_non_responders4[Z_current_non_responders4[i],],sigma_current_non_responders4[[Z_current_non_responders4[i]]],dependent.ind = 1,given.ind = c(2,3),c(D2_current[which(s==0&a1==-1)][i],D4_current[which(s==0&a1==-1)][i]),check.sigma = F)$condVar
        
        D1_cond_mean <-(a1[which(s==0&a1==-1)][i]*(beta_seq_overall[2])*(y[which(s==0&a1==-1)][i]-beta_seq_overall[1]-beta_seq_overall[3]*D2_current[which(s==0&a1==-1)][i]*a1[which(s==0&a1==-1)][i]-beta_seq_overall[5]*D4_current[which(s==0&a1==-1)][i]*a2[which(s==0&a1==-1)][i]-beta_seq_overall[6]*(X1[which(s==0&a1==-1)][i])-beta_seq_overall[7]*(X2[which(s==0&a1==-1)][i])-beta_seq_overall[8]*(X3[which(s==0&a1==-1)][i]))*cond_var[[1]]+w_non_responders*sigma_seq_overall_squared*cond_mean)/((beta_seq_overall[2])^2*cond_var+w_non_responders*sigma_seq_overall_squared)
        D1_cond_var <- w_non_responders*sigma_seq_overall_squared*cond_var/((beta_seq_overall[2])^2*cond_var+w_non_responders*sigma_seq_overall_squared)  
        
        
        if (D1_cond_var <=0) {
          D1_current[which(s==0&a1==-1)][i] <- D1_cond_mean
        } else {
          D1_current[which(s==0&a1==-1)][i] <- rtruncnorm(1,0,1,D1_cond_mean,sqrt(D1_cond_var))
        }
        
        
        
      }
    }
    X_mat_responders <- cbind(X1,X2,X3)[s==1,]
    
    #############
    #Draw pre-truncation means for each of the H mixture components for responders
    for (h in 1:H_responders) {
      if ((sum(Z_current_responders==h))>1) {
        
        
        
        eta_prop_h <- MASS::mvrnorm(1,colMeans(cbind(D1_current[which(s==1)][which(Z_current_responders==h)],
                                                     D2_current[which(s==1)][which(Z_current_responders==h)]))-colMeans(cbind(X_mat_responders[which(Z_current_responders==h),]%*%beta_D_responders1,
                                                                                                                              X_mat_responders[which(Z_current_responders==h),]%*%beta_D_responders2)),sigma_current_responders[[h]]/sum(Z_current_responders==h))
        #############
        ##############
        accept_prob <- 1    
        
        
        u <- runif(1)
        if (!is.nan(accept_prob)) {
          if (u < min(accept_prob,1,na.rm=T)){
            eta_current_responders[h,] <- eta_prop_h
          }
        }
      } else {
        
        eta_current_responders[h,] <- MASS::mvrnorm(1,c(0.5,0.5),diag(1,2))
        
      }
    }
    X_mat_non_responders3 <- cbind(X1,X2,X3)[s==0&a1==1,]
    
    
    #Draw pre-truncation means for each of the H mixture components for non-responders
    
    for (h in 1:H_non_responders) {
      if ((sum(Z_current_non_responders3==h))>1) {
        
        
        
        eta_prop_h <- MASS::mvrnorm(1,colMeans(cbind(D1_current[which(s==0&a1==1)][which(Z_current_non_responders3==h)],
                                                     D2_current[which(s==0&a1==1)][which(Z_current_non_responders3==h)],
                                                     D3_current[which(s==0&a1==1)][which(Z_current_non_responders3==h)]
        ))-colMeans(cbind(X_mat_non_responders3[which(Z_current_non_responders3==h),]%*%beta_D_non_responders31,
                          X_mat_non_responders3[which(Z_current_non_responders3==h),]%*%beta_D_non_responders32,
                          X_mat_non_responders3[which(Z_current_non_responders3==h),]%*%beta_D_non_responders33)),sigma_current_non_responders3[[h]]/sum(Z_current_non_responders3==h))
        #############
        ##############
        accept_prob <- 1    
        
        
        u <- runif(1)
        if (!is.nan(accept_prob)) {
          if (u < min(accept_prob,1,na.rm=T)){
            eta_current_non_responders3[h,] <- eta_prop_h
          }
        }
      } else {
        
        eta_current_non_responders3[h,] <- MASS::mvrnorm(1,c(0.5,0.5,0.5),diag(1,3))
        
      }
    }
    
    
    X_mat_non_responders4 <- cbind(X1,X2,X3)[s==0&a1==-1,]
    
    #Draw pre-truncation means for each of the H mixture components for non-responders
    for (h in 1:H_non_responders) {
      if ((sum(Z_current_non_responders4==h))>1) {
        
        
        
        eta_prop_h <- MASS::mvrnorm(1,colMeans(cbind(D1_current[which(s==0&a1==-1)][which(Z_current_non_responders4==h)],
                                                     D2_current[which(s==0&a1==-1)][which(Z_current_non_responders4==h)],
                                                     D4_current[which(s==0&a1==-1)][which(Z_current_non_responders4==h)]
        ))-colMeans(cbind(X_mat_non_responders4[which(Z_current_non_responders4==h),]%*%beta_D_non_responders41,
                          X_mat_non_responders4[which(Z_current_non_responders4==h),]%*%beta_D_non_responders42,
                          X_mat_non_responders4[which(Z_current_non_responders4==h),]%*%beta_D_non_responders44)),sigma_current_non_responders4[[h]]/sum(Z_current_non_responders4==h))
        #############
        ##############
        accept_prob <- 1    
        
        
        u <- runif(1)
        if (!is.nan(accept_prob)) {
          if (u < min(accept_prob,1,na.rm=T)){
            eta_current_non_responders4[h,] <- eta_prop_h
          }
        }
      } else {
        
        eta_current_non_responders4[h,] <- MASS::mvrnorm(1,c(0.5,0.5,0.5),diag(1,3))
        
      }
    }
    
    #############Compliance-covariate regression
    #########
   
    X_mat_responders <- cbind(X1,X2,X3)[s==1,][,]
    
    #Propose point
    var_temp <- crossprod((cbind(D1_current[which(s==1)],
                                 D2_current[which(s==1)])-eta_current_responders[Z_current_responders,]- matrix(cbind(X_mat_responders%*%beta_D_responders1,
                                                                                                                                   X_mat_responders%*%beta_D_responders2),nrow=length((Z_current_responders)),ncol=2,byrow=T)),((cbind(D1_current[which(s==1)],D2_current[which(s==1)])-eta_current_responders[Z_current_responders,]-matrix(cbind(X_mat_responders%*%beta_D_responders1,
                                                                                                                                                                                                                                                                                                                                                                                                        X_mat_responders%*%beta_D_responders2),nrow=length((Z_current_responders)),ncol=2,byrow=T))))
    
    sigma_prop_h <- rinvwishart(2+length(Z_current_responders), diag(1,2)+var_temp)
    
    accept_prob <- 1
    
    
    
    
    
    u <- runif(1)
    if (!is.nan(accept_prob)) {
      if (u < min(accept_prob,1,na.rm=T)) {
        sigma_current_responders1 <- sigma_prop_h
      }
    }
    
    ######
    X_mat_non_responders3 <- cbind(X1,X2,X3)[s==0&a1==1,]
    #Propose point
    var_temp <- crossprod((cbind(D1_current[which(s==0&a1==1)],
                                 D2_current[which(s==0&a1==1)],
                                 D3_current[which(s==0&a1==1)]
    )-eta_current_non_responders3[Z_current_non_responders3,]-matrix(cbind((X_mat_non_responders3)%*%beta_D_non_responders31,
                                                   (X_mat_non_responders3)%*%beta_D_non_responders32,
                                                   (X_mat_non_responders3)%*%beta_D_non_responders33),nrow=length((Z_current_non_responders3)),ncol=3,byrow=T)),((cbind(D1_current[which(s==0&a1==1)],D2_current[which(s==0&a1==1)],D3_current[which(s==0&a1==1)])-eta_current_non_responders3[Z_current_non_responders3,]-matrix( cbind((X_mat_non_responders3)%*%beta_D_non_responders31,
                                                                                                                                                                                                                                                                                                                                                                                                                                        (X_mat_non_responders3)%*%beta_D_non_responders32,
                                                                                                                                                                                                                                                                                                                                                                                                                                        (X_mat_non_responders3)%*%beta_D_non_responders33),nrow=length((Z_current_non_responders3)),ncol=3,byrow=T))))
    
    
    sigma_prop_h <- rinvwishart(3+length(Z_current_non_responders3), diag(1,3)+var_temp)
    
    accept_prob <- 1
    
    
    
    
    
    u <- runif(1)
    if (!is.nan(accept_prob)) {
      if (u < min(accept_prob,1,na.rm=T)) {
        sigma_current_non_responders1 <- sigma_prop_h
      }
    }###########
    X_mat_non_responders4 <- cbind(X1,X2,X3)[s==0&a1==-1,]
    
    #Propose point
    var_temp <- crossprod((cbind(D1_current[which(s==0&a1==-1)],
                                 D2_current[which(s==0&a1==-1)],
                                 D4_current[which(s==0&a1==-1)])-eta_current_non_responders4[Z_current_non_responders4,]-matrix(cbind((X_mat_non_responders4)%*%beta_D_non_responders41,
                                                                                                                                                   (X_mat_non_responders4)%*%beta_D_non_responders42,
                                                                                                                                                   (X_mat_non_responders4)%*%beta_D_non_responders44),nrow=length((Z_current_non_responders4)),ncol=3,byrow=T)),((cbind(D1_current[which(s==0&a1==-1)],D2_current[which(s==0&a1==-1)],D4_current[which(s==0&a1==-1)])-eta_current_non_responders4[Z_current_non_responders4,]-matrix(cbind((X_mat_non_responders4)%*%beta_D_non_responders41,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          (X_mat_non_responders4)%*%beta_D_non_responders42,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          (X_mat_non_responders4)%*%beta_D_non_responders44),nrow=length((Z_current_non_responders4)),ncol=3,byrow=T))))
    
    sigma_prop_h <- rinvwishart(3+length(Z_current_non_responders4), diag(1,3)+var_temp)
    accept_prob <- 1
    
    
    
    
    
    u <- runif(1)
    if (!is.nan(accept_prob)) {
      if (u < min(accept_prob,1,na.rm=T)) {
        sigma_current_non_responders14 <- sigma_prop_h
      }
    }
    ##############
    
    
    beta_D_responders1_mean <- solve(crossprod(X_mat_responders))%*%t(X_mat_responders) %*%( (D1_current[which(s==1)]-eta_current_responders[Z_current_responders,1])-sigma_current_responders1[1,2]*solve(sigma_current_responders1[2,2])[[1]] * ((D2_current[which(s==1)]-eta_current_responders[Z_current_responders,2])-
                                                                                                                                                                                                                                                                   (X_mat_responders%*%beta_D_responders2)))
    sigma_beta_responders_1 <- (sigma_current_responders1[1,1]-sigma_current_responders1[1,2]*solve(sigma_current_responders1[2,2])*sigma_current_responders1[1,2])[[1]]

    beta_D_responders2_mean <- solve(crossprod(X_mat_responders))%*%t(X_mat_responders) %*%( (D2_current[which(s==1)]-eta_current_responders[Z_current_responders,2])-sigma_current_responders1[1,2]*solve(sigma_current_responders1[1,1])[[1]] * ((D1_current[which(s==1)]-eta_current_responders[Z_current_responders,1])-
                                                                                                                                                                                                                                                             (X_mat_responders%*%beta_D_responders1)))
    
    sigma_beta_responders_2 <- (sigma_current_responders1[2,2]-sigma_current_responders1[1,2]*solve(sigma_current_responders1[1,1])*sigma_current_responders1[1,2])[[1]]
    
    
    
    #######
    beta_D_non_responders31_mean <- solve(crossprod(X_mat_non_responders3))%*%t(X_mat_non_responders3) %*% t((D1_current[which(s==0&a1==1)]-eta_current_non_responders3[Z_current_non_responders3,1]-sigma_current_non_responders1[1,2:3]%*%solve(sigma_current_non_responders1[2:3,2:3])%*%t(cbind(D2_current[which(s==0&a1==1)]-eta_current_non_responders3[Z_current_non_responders3,2],D3_current[which(s==0&a1==1)]-eta_current_non_responders3[Z_current_non_responders3,3])-cbind(X_mat_non_responders3%*%beta_D_non_responders32,X_mat_non_responders3%*%beta_D_non_responders33))))
    
    
    beta_D_non_responders32_mean <- solve(crossprod(X_mat_non_responders3))%*%t(X_mat_non_responders3) %*% t((D2_current[which(s==0&a1==1)]-eta_current_non_responders3[Z_current_non_responders3,2]-sigma_current_non_responders1[2,c(1,3)]%*%solve(sigma_current_non_responders1[c(1,3),c(1,3)])%*%t(cbind(D1_current[which(s==0&a1==1)]-eta_current_non_responders3[Z_current_non_responders3,1], D3_current[which(s==0&a1==1)]-eta_current_non_responders3[Z_current_non_responders3,3])-cbind(X_mat_non_responders3%*%beta_D_non_responders31,X_mat_non_responders3%*%beta_D_non_responders33))))
    beta_D_non_responders33_mean <- solve(crossprod(X_mat_non_responders3))%*%t(X_mat_non_responders3) %*% t((D3_current[which(s==0&a1==1)]-eta_current_non_responders3[Z_current_non_responders3,3]-sigma_current_non_responders1[3,1:2]%*%solve(sigma_current_non_responders1[1:2,1:2])%*%t(cbind(D1_current[which(s==0&a1==1)]-eta_current_non_responders3[Z_current_non_responders3,1],D2_current[which(s==0&a1==1)]-eta_current_non_responders3[Z_current_non_responders3,2])-cbind(X_mat_non_responders3%*%beta_D_non_responders31,X_mat_non_responders3%*%beta_D_non_responders32))))
    
    sigma_beta_non_responders_1 <- sigma_current_non_responders1[1,1]-sigma_current_non_responders1[1,2:3]%*%solve(sigma_current_non_responders1[2:3,2:3])%*%sigma_current_non_responders1[2:3,1]
    sigma_beta_non_responders_2 <- sigma_current_non_responders1[2,2]-sigma_current_non_responders1[2,c(1,3)]%*%solve(sigma_current_non_responders1[c(1,3),c(1,3)])%*%sigma_current_non_responders1[c(1,3),2]
    sigma_beta_non_responders_3 <- sigma_current_non_responders1[3,3]-sigma_current_non_responders1[3,1:2]%*%solve(sigma_current_non_responders1[1:2,1:2])%*%sigma_current_non_responders1[1:2,3]
    

    beta_D_non_responders41_mean <- solve(crossprod(X_mat_non_responders4))%*%t(X_mat_non_responders4) %*% t(D1_current[which(s==0&a1==-1)]-eta_current_non_responders4[Z_current_non_responders4,1]-sigma_current_non_responders14[1,2:3]%*%solve(sigma_current_non_responders14[2:3,2:3])%*%t(cbind(D2_current[which(s==0&a1==-1)]-eta_current_non_responders4[Z_current_non_responders4,2],D4_current[which(s==0&a1==-1)]-eta_current_non_responders4[Z_current_non_responders4,3])-cbind(X_mat_non_responders4%*%beta_D_non_responders42,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             X_mat_non_responders4%*%beta_D_non_responders44)))
    beta_D_non_responders42_mean <- solve(crossprod(X_mat_non_responders4))%*%t(X_mat_non_responders4) %*% t(D2_current[which(s==0&a1==-1)]-eta_current_non_responders4[Z_current_non_responders4,2]-sigma_current_non_responders14[2,c(1,3)]%*%solve(sigma_current_non_responders14[c(1,3),c(1,3)])%*%t(cbind(D1_current[which(s==0&a1==-1)]-eta_current_non_responders4[Z_current_non_responders4,1],D4_current[which(s==0&a1==-1)]-eta_current_non_responders4[Z_current_non_responders4,3])-cbind(X_mat_non_responders4%*%beta_D_non_responders41,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             X_mat_non_responders4%*%beta_D_non_responders44)))
    beta_D_non_responders44_mean <- solve(crossprod(X_mat_non_responders4))%*%t(X_mat_non_responders4) %*% t(D4_current[which(s==0&a1==-1)]-eta_current_non_responders4[Z_current_non_responders4,3]-sigma_current_non_responders14[3,1:2]%*%solve(sigma_current_non_responders14[1:2,1:2])%*%t(cbind(D1_current[which(s==0&a1==-1)]-eta_current_non_responders4[Z_current_non_responders4,1],D2_current[which(s==0&a1==-1)]-eta_current_non_responders4[Z_current_non_responders4,2])-cbind(X_mat_non_responders4%*%beta_D_non_responders41,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             X_mat_non_responders4%*%beta_D_non_responders42)))
    sigma_beta_non_responders_14 <- sigma_current_non_responders14[1,1]-sigma_current_non_responders14[1,2:3]%*%solve(sigma_current_non_responders14[2:3,2:3])%*%sigma_current_non_responders14[2:3,1]
    sigma_beta_non_responders_24 <- sigma_current_non_responders14[2,2]-sigma_current_non_responders14[2,c(1,3)]%*%solve(sigma_current_non_responders14[c(1,3),c(1,3)])%*%sigma_current_non_responders14[c(1,3),2]
    sigma_beta_non_responders_34 <- sigma_current_non_responders14[3,3]-sigma_current_non_responders14[3,1:2]%*%solve(sigma_current_non_responders14[1:2,1:2])%*%sigma_current_non_responders14[1:2,3]
    
    
    
    
    beta_D_responders1 <-MASS::mvrnorm(1,beta_D_responders1_mean,abs(sigma_beta_responders_1)*solve(crossprod(X_mat_responders)))
    beta_D_responders2 <- MASS::mvrnorm(1,beta_D_responders2_mean,abs(sigma_beta_responders_2)*solve(crossprod(X_mat_responders)))
    
    
    beta_D_non_responders31 <- MASS::mvrnorm(1,beta_D_non_responders31_mean,abs(sigma_beta_non_responders_1[[1]])*solve(crossprod(X_mat_non_responders3)))
    beta_D_non_responders32 <- MASS::mvrnorm(1,beta_D_non_responders32_mean,abs(sigma_beta_non_responders_2[[1]])*solve(crossprod(X_mat_non_responders3)))
    beta_D_non_responders33 <- MASS::mvrnorm(1,beta_D_non_responders33_mean,abs(sigma_beta_non_responders_3[[1]])*solve(crossprod(X_mat_non_responders3)))
    
    
    beta_D_non_responders41 <- MASS::mvrnorm(1,beta_D_non_responders41_mean,abs(sigma_beta_non_responders_14[[1]])*solve(crossprod(X_mat_non_responders4)))
    beta_D_non_responders42 <- MASS::mvrnorm(1,beta_D_non_responders42_mean,abs(sigma_beta_non_responders_24[[1]])*solve(crossprod(X_mat_non_responders4)))
    beta_D_non_responders44 <- MASS::mvrnorm(1,beta_D_non_responders44_mean,abs(sigma_beta_non_responders_34[[1]])*solve(crossprod(X_mat_non_responders4)))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #########    
    
    #####
    
    #Draw pre-truncation covariance matrices for each of the H mixture components for responders
    for (h in 1:H_responders) {
      if (sum(Z_current_responders==h)>1) {
        X_mat_responders <- cbind(X1,X2,X3)[s==1,][which(Z_current_responders==h),]
        
        #Propose point
        var_temp <- crossprod((cbind(D1_current[which(s==1)][which(Z_current_responders==h)],
                                     D2_current[which(s==1)][which(Z_current_responders==h)])-matrix(eta_current_responders[h,],nrow=length(which(Z_current_responders==h)),ncol=2,byrow=T)- matrix(cbind(X_mat_responders%*%beta_D_responders1,
                                                                                                                                       X_mat_responders%*%beta_D_responders2),nrow=length(which(Z_current_responders==h)),ncol=2,byrow=T)),((cbind(D1_current[which(s==1)][which(Z_current_responders==h)],D2_current[which(s==1)][which(Z_current_responders==h)])-matrix(eta_current_responders[h,],nrow=length(which(Z_current_responders==h)),ncol=2,byrow=T)-matrix(cbind(X_mat_responders%*%beta_D_responders1,
                                                                                                                                                                                                                                                                                                                                                                                                            X_mat_responders%*%beta_D_responders2),nrow=length(which(Z_current_responders==h)),ncol=2,byrow=T))))
        
        sigma_prop_h <- rinvwishart(2+sum(Z_current_responders==h), diag(1,2)+var_temp)
        
        accept_prob <- 1
        
        
        
        
        
        u <- runif(1)
        if (!is.nan(accept_prob)) {
          if (u < min(accept_prob,1,na.rm=T)) {
            sigma_current_responders[[h]] <- sigma_prop_h
          }
        }
      }
      else {
        sigma_current_responders[[h]] <- (rinvwishart(2,diag(1,2)))  
      }
    } 
    
    
    
    #Draw pre-truncation covariance matrices for each of the H mixture components for non-responders
    
    for (h in 1:H_non_responders) {
      if (sum(Z_current_non_responders3==h)>1) {
        X_mat_non_responders3 <- cbind(X1,X2,X3)[s==0&a1==1,][which(Z_current_non_responders3==h),]
        #Propose point
        var_temp <- crossprod((cbind(D1_current[which(s==0&a1==1)][which(Z_current_non_responders3==h)],
                                     D2_current[which(s==0&a1==1)][which(Z_current_non_responders3==h)],
                                     D3_current[which(s==0&a1==1)][which(Z_current_non_responders3==h)]
        )-matrix(eta_current_non_responders3[h,],nrow=length(which(Z_current_non_responders3==h)),ncol=3,byrow=T)-matrix(cbind((X_mat_non_responders3)%*%beta_D_non_responders31,
                                                       (X_mat_non_responders3)%*%beta_D_non_responders32,
                                                       (X_mat_non_responders3)%*%beta_D_non_responders33),nrow=length(which(Z_current_non_responders3==h)),ncol=3,byrow=T)),((cbind(D1_current[which(s==0&a1==1)][which(Z_current_non_responders3==h)],D2_current[which(s==0&a1==1)][which(Z_current_non_responders3==h)],D3_current[which(s==0&a1==1)][which(Z_current_non_responders3==h)])-matrix(eta_current_non_responders3[h,],nrow=length(which(Z_current_non_responders3==h)),ncol=3,byrow=T)-matrix( cbind((X_mat_non_responders3)%*%beta_D_non_responders31,
                                                                                                                                                                                                                                                                                                                                                                                                                                            (X_mat_non_responders3)%*%beta_D_non_responders32,
                                                                                                                                                                                                                                                                                                                                                                                                                                            (X_mat_non_responders3)%*%beta_D_non_responders33),nrow=length(which(Z_current_non_responders3==h)),ncol=3,byrow=T))))
        
        
        sigma_prop_h <- rinvwishart(3+sum(Z_current_non_responders3==h), diag(1,3)+var_temp)
        
        accept_prob <- 1
        
        
        
        
        
        u <- runif(1)
        if (!is.nan(accept_prob)) {
          if (u < min(accept_prob,1,na.rm=T)) {
            sigma_current_non_responders3[[h]] <- sigma_prop_h
          }
        }
      }
      else {
        sigma_current_non_responders3[[h]] <- (rinvwishart(3,diag(1,3)))  
      }
    } 
    
    #Draw pre-truncation covariance matrices for each of the H mixture components for non-responders
    
    for (h in 1:H_non_responders) {
      if (sum(Z_current_non_responders4==h)>1) {
        X_mat_non_responders4 <- cbind(X1,X2,X3)[s==0&a1==-1,][which(Z_current_non_responders4==h),]
        
        #Propose point
        var_temp <- crossprod((cbind(D1_current[which(s==0&a1==-1)][which(Z_current_non_responders4==h)],
                                     D2_current[which(s==0&a1==-1)][which(Z_current_non_responders4==h)],
                                     D4_current[which(s==0&a1==-1)][which(Z_current_non_responders4==h)])-matrix(eta_current_non_responders4[h,],nrow=length(which(Z_current_non_responders4==h)),ncol=3,byrow=T)-matrix(cbind((X_mat_non_responders4)%*%beta_D_non_responders41,
                                                                                                                                                       (X_mat_non_responders4)%*%beta_D_non_responders42,
                                                                                                                                                       (X_mat_non_responders4)%*%beta_D_non_responders44),nrow=length(which(Z_current_non_responders4==h)),ncol=3,byrow=T)),((cbind(D1_current[which(s==0&a1==-1)][which(Z_current_non_responders4==h)],D2_current[which(s==0&a1==-1)][which(Z_current_non_responders4==h)],D4_current[which(s==0&a1==-1)][which(Z_current_non_responders4==h)])-matrix(eta_current_non_responders4[h,],nrow=length(which(Z_current_non_responders4==h)),ncol=3,byrow=T)-matrix(cbind((X_mat_non_responders4)%*%beta_D_non_responders41,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              (X_mat_non_responders4)%*%beta_D_non_responders42,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              (X_mat_non_responders4)%*%beta_D_non_responders44),nrow=length(which(Z_current_non_responders4==h)),ncol=3,byrow=T))))
        
        sigma_prop_h <- rinvwishart(3+sum(Z_current_non_responders4==h), diag(1,3)+var_temp)
        accept_prob <- 1
        
        
        
        
        
        u <- runif(1)
        if (!is.nan(accept_prob)) {
          if (u < min(accept_prob,1,na.rm=T)) {
            sigma_current_non_responders4[[h]] <- sigma_prop_h
          }
        }
      }
      else {
        sigma_current_non_responders4[[h]] <- (rinvwishart(3,diag(1,3)))
      }
    } 
    
    ###############################
    
    #Draw w primes
    W_prime_responders <- computeWprime(Z_current_responders,H_responders,alpha_current_responders)
    W_prime_non_responders3 <- computeWprime(Z_current_non_responders3,H_non_responders,alpha_current_non_responders3)
    W_prime_non_responders4 <- computeWprime(Z_current_non_responders4,H_non_responders,alpha_current_non_responders4)
    
    #Draw mixture weights
    W_current_responders <- computeW(W_prime_responders,H_responders)
    W_current_non_responders3<- computeW(W_prime_non_responders3,H_non_responders)
    W_current_non_responders4<- computeW(W_prime_non_responders4,H_non_responders)
    
    #Draw concentration parameter
    alpha_current_responders <- computeAlpha(W_prime_responders,Z_current_responders,alpha_current_responders,H_responders)
    alpha_current_non_responders3 <- computeAlpha(W_prime_non_responders3,Z_current_non_responders3,alpha_current_non_responders3,H_non_responders)
    alpha_current_non_responders4 <- computeAlpha(W_prime_non_responders4,Z_current_non_responders4,alpha_current_non_responders4,H_non_responders)
      
    
    X_mat_non_responders3 <- cbind(X1,X2,X3)[s==0&a1==1,]
    X_mat_non_responders4 <- cbind(X1,X2,X3)[s==0&a1==-1,]
    X_mat_responders <- cbind(X1,X2,X3)[s==1,]
    Z_current_non_responders3 <- computeZnonresponders(length(which(s==0&a1==1)),H_non_responders,W_current_non_responders3,D1_current[which(s==0&a1==1)],
                                                       D2_current[which(s==0&a1==1)],
                                                       D3_current[which(s==0&a1==1)],
                                                       eta_current_non_responders3,
                                                       beta_D_non_responders31,beta_D_non_responders32,beta_D_non_responders33,X_mat_non_responders3,sigma_current_non_responders3)
    
    Z_current_non_responders4 <- computeZnonresponders(length(which(s==0&a1==-1)),H_non_responders,W_current_non_responders4,D1_current[which(s==0&a1==-1)],
                                                       D2_current[which(s==0&a1==-1)],
                                                       D4_current[which(s==0&a1==-1)],
                                                       eta_current_non_responders4,
                                                       beta_D_non_responders41,beta_D_non_responders42,beta_D_non_responders44,X_mat_non_responders4,sigma_current_non_responders4)
    Z_current_responders <- computeZresponders(length(which(s==1)),H_responders,W_current_responders,D1_current[which(s==1)],D2_current[which(s==1)],eta_current_responders,beta_D_responders1,beta_D_responders2,X_mat_responders,sigma_current_responders)
    
    ##########
    
    
    ##################The following is fitting the parametric outcome marginal structural model.
    #######
    
    #a2[s==1] <- 0
    X_seq_overall <- data.frame(cbind(1,D1_current,D2_current,D3_current,D4_current,X1,X2,X3,D1_current,D2_current,D3_current,D4_current,a1,a2,s,y))
    
    X_seq_overall <- X_seq_overall %>% filter(X_seq_overall$s == 1) %>% rbind(X_seq_overall %>% filter(X_seq_overall$s == 1) %>% mutate(a2=-a2), X_seq_overall %>% filter(X_seq_overall$s ==0))
    
    X_seq_overall[,c(4)] <- (X_seq_overall[,13]==1)*X_seq_overall[,14]*X_seq_overall[,c(4)]
    X_seq_overall[,c(5)] <- (X_seq_overall[,13]==-1)*X_seq_overall[,14]*X_seq_overall[,c(5)]
    X_seq_overall[,2] <- X_seq_overall[,2] * X_seq_overall[,13]
    X_seq_overall[,3] <- X_seq_overall[,3] * X_seq_overall[,13]
    y2 <- X_seq_overall[,16]
    s2 <- X_seq_overall$s
    w_residual <- (X_seq_overall$s==1)/0.5+
      (X_seq_overall$s==0)/0.25
    
    
    w_residual2 <- (X_seq_overall$s==1)/(0.5) + 
      (X_seq_overall$s==0)/(0.25)
    
    X_seq_overall <- as.matrix(X_seq_overall[,c(1:8)])
    
    
    
    
    
    beta_seq_mean_overall <- coef(lm(y2~X_seq_overall-1,weights = w_residual2))

    #Compute RSS for drawing from outcome variance
    RSS_seq_1 <- t(y2-X_seq_overall%*%beta_seq_overall)%*%diag(w_residual2)%*%(y2-X_seq_overall%*%beta_seq_overall)
    
    if(!is.na(RSS_seq_1)) {
      sigma_seq_overall_squared <- rinvchisq(1,length(y2)-8,RSS_seq_1/(length(y2)-8))
      
    }
    
    
    beta_seq_var_overall <- sigma_seq_overall_squared*solve(t(X_seq_overall)%*%diag(w_residual2)%*%X_seq_overall)
    
    
    
    beta_seq_overall <- mvrnorm(1,beta_seq_mean_overall,beta_seq_var_overall)
    
    
    beta_seq_overall_results[j,] <-t(beta_seq_overall)
    
    
    #####
    
    D1_matrix[j,] <- D1_current
    D2_matrix[j,] <- D2_current
    D3_matrix[j,] <- D3_current
    D4_matrix[j,] <- D4_current
    
    sigma_seq_1_squared_results[j] <- sigma_seq_overall_squared
    
    print(table(Z_current_responders))
    print(table(Z_current_non_responders3))
    print(table(Z_current_non_responders4))
    print(W_current_responders)
    print(W_current_non_responders3)
    print(beta_seq_overall)
    print(j)
    print(sigma_seq_overall_squared)
  }
  
  
  return(list(beta_seq_overall_results, #Regression coefficients for marginal structural model
              sigma_seq_1_squared_results#, #Residual variance estimate

  ))
}

