
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
    GH_baseline <- X[,1]
    higheduce <- X[,2]
    smoke_baseline <- X[,3]
    ####
    
    #Simulate potential compliances from a Gaussian copula. D11, D12, D21, D22 are potential compliances
    cop <- normalCopula(param=c(0.2), dim = 2, dispstr = "un")
    
    D <- matrix(NA,nrow=n,ncol=2)
    
    
    for (j in 1:n) {
      mvd <- mvdc(copula=cop, margins=c("truncnorm","truncnorm"),
                  paramMargins = list(list(mean=(0.5*GH_baseline[j]+0.5*smoke_baseline[j]),sd=0.5,a=0,b=1),
                                      list(mean=(0.5*higheduce[j]),sd=0.5,a=0,b=1)))
      
      #Draw from Gaussian Copula
      D[j,] <- rMvdc(1,mvd)
    }
    
    #Stage-1 potential compliances
    D11 <- D[,1] #Stage 1: A1=+1
    D12 <- D[,2] #Stage 1: A1=-1 
    
    
    #First stage response indicator (also referred to as R in manuscript)
    s <- rep(NA,n)
    
    s <- rbinom(n,1,expit(a1*(0.5*D11+0.7*D12)-0.5+0.2*GH_baseline-0.2*smoke_baseline))
    
    
    
    
    D <- matrix(NA,nrow=n,ncol=2)
    for (j in 1:n) {
      mvd <- mvdc(copula=cop, margins=c("truncnorm","truncnorm"),
                  paramMargins = list(list(mean=(0.2*D11[j]+0.3*D12[j]-0.2*GH_baseline[j]+0.5*smoke_baseline[j]),sd=0.5,a=0,b=1),
                                      list(mean=(0.5*D11[j]+0.4*D12[j]+0.5*higheduce[j]),sd=0.5,a=0,b=1)))
      
      #Draw from Gaussian Copula
      D[j,] <- rMvdc(1,mvd)
    }
    
    D21 <- (1-s)*D[,1] #Stage 2: A2 =+1
    D22 <- (1-s)*D[,2]
    
    
    
    
    ############# Create treatment sequence-wise outcomes without noise
    y1 <- rep(NA,n)
    
    
    
    y1 <- 2.5 +a1*(-0.8*D11+0.4*D12)+a2*(-0.2*D21*(a1!=-1)+0.5*D22*(a1!=1))*(1-s)+0.3*GH_baseline+0.6*higheduce-0.7*smoke_baseline
    
    #y <- exp(y1 + rnorm(n,0,.1) )
    y <- MASS::rnegbin(n,exp(y1),0.4)#rnegbin(n,100,10)#rpois(n,10+exp(y1))# + rnorm(n,0,.1)
    #y <- y*rbinom(n,1,0.7)
    y <- log(y+0.5)
    
    #y <- log(ifelse(y==0,2,y))
    sim_list[[i]] <- data.frame(a1,a2,s,D11,D12,D21,D22,y,y1,GH_baseline,higheduce,smoke_baseline)
    
    print(i)
  }
  sim_list
}



set.seed(136456)
dat_for_analysis_4 <- ComplianceSimulate(n=160,n_sim=1)[[1]]


