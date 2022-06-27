MSM_200_without_covariates <- c(MSM_200_10_17_21_H_5_3_5000_without_covariates_beta_1_200,
                                MSM_200_10_17_21_H_5_3_5000_without_covariates_beta_201_400)


MSM_400_without_covariates <- c(MSM_400_10_17_21_H_5_5_5000_without_covariates_beta_1_200,
                             MSM_400_10_17_21_H_5_5_5000_without_covariates_beta_1_201_400)

library(tidyverse)



plot(1:5000,MSM_400_without_covariates[[16]][[1]][,2],type="l")


###########
mean_SD_400_beta <- abs(round(apply(do.call(rbind,lapply((1:400),function(x) (apply(c(MSM_400_without_covariates)[[x]][[1]][1:5000,],2,mad)))),2,mean),2))

MC_SE_400_beta <- abs(round(apply(do.call(rbind,lapply((1:400),function(x) (apply(c(MSM_400_without_covariates)[[x]][[1]][1:5000,],2,median)))),2,mad),2))

bias_400_beta <- abs(round(colMeans(do.call(rbind,lapply(1:400,function(x) (apply(MSM_400_without_covariates[[x]][[1]][,],2,median)-c(0.5,-0.8,0.4,-0.2,0.5,0.3,0.6,-0.7))))),2))

#############
mean_SD_200_beta <- abs(round(apply(do.call(rbind,lapply(setdiff((1:400),c(4,114)),function(x) (apply(c(MSM_200_without_covariates)[[x]][[1]][1:5000,],2,mad)))),2,mean),2))

MC_SE_200_beta <- abs(round(apply(do.call(rbind,lapply(setdiff((1:400),c(4,114)),function(x) (apply(c(MSM_200_without_covariates)[[x]][[1]][1:5000,],2,median)))),2,mad),2))

bias_200_beta <- abs(round(colMeans(do.call(rbind,lapply(setdiff((1:400),c(4,114)),function(x) (apply(MSM_200_without_covariates[[x]][[1]][,],2,median)-c(0.5,-0.8,0.4,-0.2,0.5,0.3,0.6,-0.7))))),2))
##########

beta_sumamry_table <- matrix(NA,nrow=8,ncol=6)


beta_sumamry_table <- cbind(bias_200_beta,
      MC_SE_200_beta,
      mean_SD_200_beta,
      bias_400_beta,
      MC_SE_400_beta,
      mean_SD_400_beta)


colnames(beta_sumamry_table) <- c("Bias 200",
                                       "MCSE 200",
                                       "SD 200",
                                       "Bias 400",
                                       "MCSE 400",
                                       "SD 400")




knitr::kable(cbind(beta_sumamry_table),format="latex")
