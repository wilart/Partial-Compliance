# Generates a list of simulated datasets with n people and n_sim simulated datasets
# for an ENGAGE-type study where the compliance distribution for stage-1 potential compliances is generated from a copula with beta distributed margins 
# The stage-2 potential compliances are a transformation of potential compliances generated from a copula model with beta margins such that
# it is a function of stage-1 potential compliances.
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
    
    X <- MASS::mvrnorm(n,c(-0.5,0,0.5),diag(c(0.8,0.8,0.8),3)%*%rbind(c(1,0.2,0.2),
                                                                     c(0.2,1,0.2),
                                                                     c(0.2,0.2,1))%*%diag(c(0.8,0.8,0.8),3))
    X1 <- X[,1]
    X2 <- X[,2]
    X3 <- X[,3]
    
    #Simulate potential compliances from a Gaussian copula. D1, D2, D3, D4 are the potential compliances
    cop <- normalCopula(param=c(0.3), dim = 2, dispstr = "un")
    
    D <- matrix(NA,nrow=n,ncol=2)
    for (j in 1:n) {
      mvd <- mvdc(copula=cop, margins=c("beta","beta"),
                  paramMargins = list(list(shape1=0.5,shape2=1),
                                      list(shape1=1,shape2=2)))
      
      #Draw from Gaussian Copula
      D[j,] <- rMvdc(1,mvd)
    }
    
    D1 <- D[,1] #Stage 1: A1=+1
    D2 <- D[,2] #Stage 1: A1=-1 
    
    #First stage response indicators
    s <- rep(NA,n)
    
    s <- rbinom(n,1,expit(a1*(0.5*D1+0.7*D2)-0.5+0.2*X1-0.2*X3))

    
    
    
    D <- matrix(NA,nrow=n,ncol=2)
    for (j in 1:n) {
      mvd <- mvdc(copula=cop, margins=c("beta","beta"),
                  paramMargins = list(list(shape1=0.5,shape2=0.5),
                                      list(shape=0.7,1)))
      
      #Draw from Gaussian Copula
      D[j,] <- rMvdc(1,mvd)
    }
    
    D3 <- (1-s)*D[,1]*expit(D1+D2) 
    D4 <- (1-s)*D[,2]*expit(D1-D2) 
    
    
   
    ############# Create treatment sequence-wise outcomes without noise
    y1 <- rep(NA,n)
    
   
    
    
    
    y1 <- 0.5 +a1*(-0.8*D1+0.4*D2)+a2*(-0.2*D3*(a1!=-1)+0.5*D4*(a1!=1))*(1-s)+0.3*X1+0.6*X2-0.7*X3
    
    #Final outcome with noise.
    y <- y1 + rnorm(n,0,0.1) 
    
    
    sim_list[[i]] <- data.frame(a1,a2,s,D1,D2,D3,D4,y,y1,X1,X2,X3)
    
    print(i)
  }
  sim_list
}

set.seed(38265)
sim_200 <- ComplianceSimulate(n=200, n_sim=400)

set.seed(382651)
sim_400 <- ComplianceSimulate(n=400, n_sim=400)
