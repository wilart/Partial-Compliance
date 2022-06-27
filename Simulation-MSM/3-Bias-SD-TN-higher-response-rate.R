MSM_200_with_covariates <- c(MSM_200_10_27_21_H_5_3_5000_with_covariates_beta_higher_response_rate_1_200,
                             MSM_200_10_27_21_H_5_3_5000_with_covariates_beta_higher_response_rate_201_400)



library(tidyverse)



plot(1:5000,MSM_200_with_covariates[[16]][[1]][,2],type="l")


###########
mean_SD_200_truncnorm <- abs(round(apply(do.call(rbind,lapply(1:400,function(x) (apply(c(MSM_200_with_covariates)[[x]][[1]][1:5000,],2,mad)))),2,mean),2))

MC_SE_200_truncnorm <- abs(round(apply(do.call(rbind,lapply(1:400,function(x) (apply(c(MSM_200_with_covariates)[[x]][[1]][1:5000,],2,median)))),2,mad),2))

bias_200_truncnorm <- abs(round(colMeans(do.call(rbind,lapply(1:400,function(x) (apply(MSM_200_with_covariates[[x]][[1]][,],2,median)-c(0.5,-0.8,0.4,-0.2,0.5,0.3,0.6,-0.7))))),2))

#############

truncnorm_sumamry_table <- matrix(NA,nrow=8,ncol=6)


truncnorm_sumamry_table <- rbind(bias_200_truncnorm,
                                 MC_SE_200_truncnorm,
                                 mean_SD_200_truncnorm)


rownames(truncnorm_sumamry_table) <- c("Bias 200",
                                       "MCSE 200",
                                       "SD 200")




knitr::kable(cbind(truncnorm_sumamry_table),format="latex")
