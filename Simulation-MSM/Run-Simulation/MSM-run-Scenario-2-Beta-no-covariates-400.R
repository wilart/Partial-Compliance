#Fit model



library(doParallel)
library(doRNG)
expit <- function(x) exp(x)/(1+exp(x))

#Number of cores/threads
no_cores <- detectCores()-2

cl<-makeCluster(no_cores)
registerDoParallel((cl))
Sys.time()


tst_0.5_interaction_start <- Sys.time()

set.seed(2341)

MSM_400_10_17_21_H_5_5_5000_without_covariates_beta_1_200<- foreach(i =1:200,.errorhandling = "pass",.packages=c("dplyr","tmvtnorm",
                                                                                                                "truncnorm",
                                                                                                                "LaplacesDemon",
                                                                                                                "MASS",
                                                                                                                "condMVNorm","tidyverse")) %dorng% computeBetas(5000,sim_400[[i]],H_responders =5,H_non_responders =5)
stopCluster(cl)
Sys.time()
tst_0.5_interaction_stop <- Sys.time()



#Number of cores/threads
no_cores <- detectCores()-2

cl<-makeCluster(no_cores)
registerDoParallel((cl))
Sys.time()


tst_0.5_interaction_start <- Sys.time()

set.seed(23452)

MSM_400_10_17_21_H_5_5_5000_without_covariates_beta_1_201_400<- foreach(i =201:400,.errorhandling = "pass",.packages=c("dplyr","tmvtnorm",
                                                                                                                      "truncnorm",
                                                                                                                      "LaplacesDemon",
                                                                                                                      "MASS",
                                                                                                                      "condMVNorm","tidyverse")) %dorng% computeBetas(5000,sim_400[[i]],H_responders =5,H_non_responders =5)
stopCluster(cl)
Sys.time()
tst_0.5_interaction_stop <- Sys.time()
