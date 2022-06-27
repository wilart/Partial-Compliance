

MSM_160_zero_inflated <-c(MSM_160_10_28_21_H_3_3_5000_with_covariates_beta_zero_inflated_1_200,
                          MSM_160_10_28_21_H_3_3_5000_with_covariates_beta_zero_inflated_1_201_400)
library(tidyverse)


###########
#############
mean_SD_160_truncnorm <- abs(round(apply(do.call(rbind,lapply(1:400,function(x) (apply(c(MSM_160_zero_inflated)[[x]][[1]][1:5000,],2,mad)))),2,mean),2))

MC_SE_160_truncnorm <- abs(round(apply(do.call(rbind,lapply(1:400,function(x) (apply(c(MSM_160_zero_inflated)[[x]][[1]][1:5000,],2,median)))),2,mad),2))

bias_160_truncnorm <- abs(round(colMeans(do.call(rbind,lapply(1:400,function(x) (apply(MSM_160_zero_inflated[[x]][[1]][,],2,median)-zero_inflated_median)))),2))

##########

truncnorm_summary_table <- matrix(NA,nrow=8,ncol=3)


truncnorm_summary_table <- rbind(bias_160_truncnorm,
                                 MC_SE_160_truncnorm,
                                 mean_SD_160_truncnorm)


rownames(truncnorm_summary_table) <- c("Bias 160",
                                       "MCSE 160",
                                       "SD 160")




knitr::kable(cbind(truncnorm_summary_table),format="latex")
