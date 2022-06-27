
computeDifference <- function(thetadraws) {
  upper_limit <- rep(NA,4)

  max_ind <- which.min(colMeans(thetadraws))
  
  diff_matrix <- thetadraws[,]-matrix(thetadraws[,max_ind],nrow=nrow(thetadraws),ncol=4)
  return(diff_matrix)
}



MSM_200_with_covariates_1_400 <- c(MSM_200_10_25_21_H_3_3_5000_with_covariates_beta_1_200,
                                   MSM_200_10_25_21_H_3_3_5000_with_covariates_beta_1_201_400)
########



library(doParallel)

#Number of cores/threads
no_cores <- detectCores()-2

cl<-makeCluster(no_cores)
registerDoParallel((cl))
Sys.time()

set.seed(1274)
clusterSetRNGStream(cl, 123)

tst_0.5_interaction_start <- Sys.time()

dat_200_median <- foreach(x =1:400,.errorhandling = "pass",.packages=c("dplyr","tmvtnorm",
                                                                                "truncnorm",
                                                                                "LaplacesDemon",
                                                                                "MASS",
                                                                                "condMVNorm","tidyverse")) %dopar% (cbind(apply(MSM_200_with_covariates_1_400[[x]][[1]][1:5000,],1,function(z) z[1]+z[2]*0.5+z[3]*0.5+z[4]*0.5+z[6]*median(sim_1truncnorm_200[[x]]$X1[sim_1truncnorm_200[[x]]$a1==1&sim_1truncnorm_200[[x]]$a2==1])+z[7]*median(sim_1truncnorm_200[[x]]$X2[sim_1truncnorm_200[[x]]$a1==1&sim_1truncnorm_200[[x]]$a2==1])+z[8]*median(sim_1truncnorm_200[[x]]$X3[sim_1truncnorm_200[[x]]$a1==1&sim_1truncnorm_200[[x]]$a2==1])),
                                                                                                                          apply(MSM_200_with_covariates_1_400[[x]][[1]][1:5000,],1,function(z) z[1]+z[2]*0.5+z[3]*0.5-z[4]*0.5+z[6]*median(sim_1truncnorm_200[[x]]$X1[sim_1truncnorm_200[[x]]$a1==1&sim_1truncnorm_200[[x]]$a2==-1])+z[7]*median(sim_1truncnorm_200[[x]]$X2[sim_1truncnorm_200[[x]]$a1==1&sim_1truncnorm_200[[x]]$a2==-1])+z[8]*median(sim_1truncnorm_200[[x]]$X3[sim_1truncnorm_200[[x]]$a1==1&sim_1truncnorm_200[[x]]$a2==-1])),
                                                                                                                          apply(MSM_200_with_covariates_1_400[[x]][[1]][1:5000,],1,function(z) z[1]-z[2]*0.5-z[3]*0.5+z[5]*0.5+z[6]*median(sim_1truncnorm_200[[x]]$X1[sim_1truncnorm_200[[x]]$a1==-1&sim_1truncnorm_200[[x]]$a2==1])+z[7]*median(sim_1truncnorm_200[[x]]$X2[sim_1truncnorm_200[[x]]$a1==-1&sim_1truncnorm_200[[x]]$a2==1])+z[8]*median(sim_1truncnorm_200[[x]]$X3[sim_1truncnorm_200[[x]]$a1==-1&sim_1truncnorm_200[[x]]$a2==1])),
                                                                                                                          apply(MSM_200_with_covariates_1_400[[x]][[1]][1:5000,],1,function(z) z[1]-z[2]*0.5-z[3]*0.5-z[5]*0.5+z[6]*median(sim_1truncnorm_200[[x]]$X1[sim_1truncnorm_200[[x]]$a1==-1&sim_1truncnorm_200[[x]]$a2==-1])+z[7]*median(sim_1truncnorm_200[[x]]$X2[sim_1truncnorm_200[[x]]$a1==-1&sim_1truncnorm_200[[x]]$a2==-1])+z[8]*median(sim_1truncnorm_200[[x]]$X3[sim_1truncnorm_200[[x]]$a1==-1&sim_1truncnorm_200[[x]]$a2==-1]))))

stopCluster(cl)
Sys.time()
tst_0.5_interaction_stop <- Sys.time()
######################
computeUpperLimitsUnadjustedMCB <- function(thetadraws) {
  upper_limit <- rep(NA,4)
  
  #Compute index of best EDTR
  max_odds_ind <- which.min(colMeans(thetadraws))
  
  #Compute differences between each EDTR and best
  diff_matrix <- matrix(thetadraws[,max_odds_ind],nrow=nrow(thetadraws),ncol=4)-thetadraws
  return(apply(diff_matrix,2,quantile,1-0.05))
}
upper_limits <- lapply(dat_200_median,function(x)ComputeMCBUpperLimits(x))
upper_limits_2 <- lapply(dat_200_median,function(x)computeUpperLimitsUnadjustedMCB(x))

#Bayesian
c("Bayesian", round(mean(unlist(lapply(1:400,function(x) mean(apply(matrix((upper_limits[[x]]),nrow=5000,ncol=4,byrow=T)>=-computeDifference((dat_200_median[[x]][seq(1,5000,1),])),1,prod))))),4))

c("Bayesian", round(mean(unlist(lapply(1:400,function(x) mean(apply(matrix((upper_limits_2[[x]]),nrow=5000,ncol=4,byrow=T)>=-computeDifference((dat_200_median[[x]][seq(1,5000,1),])),1,prod))))),4))
