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

set.seed(123)

MSM_200_10_27_21_H_5_3_5000_with_covariates_beta_higher_response_rate_1_200<- foreach(i =1:200,.errorhandling = "pass",.packages=c("dplyr","tmvtnorm",
                                                                                       "truncnorm",
                                                                                       "LaplacesDemon",
                                                                                       "MASS",
                                                                                       "condMVNorm","tidyverse")) %dorng% computeBetas(5000,sim_1truncnorm_200_higher_response_rate[[i]],H_responders =5,H_non_responders =3)
stopCluster(cl)
Sys.time()
tst_0.5_interaction_stop <- Sys.time()

#Number of cores/threads
no_cores <- detectCores()-2

cl<-makeCluster(no_cores)
registerDoParallel((cl))
Sys.time()


tst_0.5_interaction_start <- Sys.time()

set.seed(1234)

MSM_200_10_27_21_H_5_3_5000_with_covariates_beta_higher_response_rate_201_400<- foreach(i =201:400,.errorhandling = "pass",.packages=c("dplyr","tmvtnorm",
                                                                                                             "truncnorm",
                                                                                                             "LaplacesDemon",
                                                                                                             "MASS",
                                                                                                             "condMVNorm","tidyverse")) %dorng% computeBetas(5000,sim_1truncnorm_200_higher_response_rate[[i]],H_responders =5,H_non_responders =3)
stopCluster(cl)
Sys.time()
tst_0.5_interaction_stop <- Sys.time()

