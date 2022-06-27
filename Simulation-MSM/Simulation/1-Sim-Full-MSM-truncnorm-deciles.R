# Generates a list of simulated datasets with n people and n_sim simulated datasets
# for an ENGAGE-type study such that the stage-1 potential compliances are generated from a Gaussian copula with truncated normal margins and the observed compliances are deciles
# which are functions of the baseline covariates.
library(copula)
library(truncnorm)

ComplianceSimulate <- function(n, n_sim){
  
  sim_list <- vector("list", length = n_sim)
  
  for (i in 1:n_sim){
    
    a1 <- 2*rbinom(n,1,.5)-1 #First stage treatment
        a2 <- 2*rbinom(n,1,0.5)-1 #Second stage treatment
    
    
    expit <- function(x)exp(x)/(1+exp(x))
    
    ####Baseline covariates
    X <- matrix(NA,nrow=n,ncol=3)
    
    X<-MASS::mvrnorm(n,c(-0.5,0,0.5),diag(c(0.8,0.8,0.8),3)%*%rbind(c(1,0.2,0.2),
                                                                    c(0.2,1,0.2),
                                                                    c(0.2,0.2,1))%*%diag(c(0.8,0.8,0.8),3))
    X1 <- X[,1]
    X2 <- X[,2]
    X3 <- X[,3]
    ####
    
    #Simulate potential compliances from a Gaussian copula. D1, D2, D3, D4 are potential compliances
    cop <- normalCopula(param=c(0.2), dim = 2, dispstr = "un")
    
    D <- matrix(NA,nrow=n,ncol=2)
    
    
    for (j in 1:n) {
      mvd <- mvdc(copula=cop, margins=c("truncnorm","truncnorm"),
                  paramMargins = list(list(mean=(0.5*X1[j]+0.5*X3[j]),sd=0.5,a=0,b=1),
                                      list(mean=(0.5*X2[j]),sd=0.5,a=0,b=1)))
      
      #Draw from Gaussian Copula
      D[j,] <- rMvdc(1,mvd)
    }
    
    #Stage-1 potential compliances
    D1 <- D[,1] #Stage 1: A1=+1
    D2 <- D[,2] #Stage 1: A1=-1 
    
    
    #First stage response indicator (also referred to as R in manuscript)
    s <- rep(NA,n)
    
    s <- rbinom(n,1,expit(a1*(0.5*D1+0.7*D2)-0.5+0.2*X1-0.2*X3))

    
    
    
    D <- matrix(NA,nrow=n,ncol=2)
    for (j in 1:n) {
      mvd <- mvdc(copula=cop, margins=c("truncnorm","truncnorm"),
                  paramMargins = list(list(mean=(0.2*D1[j]+0.3*D2[j]-0.2*X1[j]+0.5*X3[j]),sd=0.5,a=0,b=1),
                                      list(mean=(0.5*D1[j]+0.4*D2[j]+0.5*X2[j]),sd=0.5,a=0,b=1)))
      
      #Draw from Gaussian Copula
      D[j,] <- rMvdc(1,mvd)
    }
    
    D3 <- (1-s)*D[,1] #Stage 2: A2 =+1
    D4 <- (1-s)*D[,2]
    
    
    
    
    ############# Create treatment sequence-wise outcomes without noise
    y1 <- rep(NA,n)
    

    
    y1 <- 0.5 +a1*(-0.8*D1+0.4*D2)+a2*(-0.2*D3*(a1!=-1)+0.5*D4*(a1!=1))*(1-s)+0.3*X1+0.6*X2-0.7*X3
   
    D1 <- round(D1,1)
    D2 <- round(D2,1)
    D3 <- round(D3,1)
    D4 <- round(D4,1)
    
    y <- y1 + rnorm(n,0,.1) 
    sim_list[[i]] <- data.frame(a1,a2,s,D1,D2,D3,D4,y,y1,X1,X2,X3)
    
    print(i)
  }
  sim_list
}

set.seed(382651)
sim_1truncnorm_200_deciles<- ComplianceSimulate(n=200,n_sim=400)
