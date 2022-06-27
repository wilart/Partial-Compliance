complianceMatrix <- function(D1_plus,D1_minus,D2,dat,var_dat) {
  
  outcome_MSM_1 <- outcome_MSM_2 <- outcome_MSM_3 <- outcome_MSM_4 <- rep(NA,10000)
  
  outcome_MSM_1 <- (dat[,]%*%c(1,
                               D1_plus,
                               D1_minus,
                               D2,
                               0*D2,
                               mean(dat_for_analysis_4$GH_baseline),
                               mean(dat_for_analysis_4$higheduce),
                               mean(dat_for_analysis_4$smoke_baseline)))
  
  
  
  outcome_MSM_2 <- (dat[,]%*%  c(1,
                                 D1_plus,
                                 D1_minus,
                                 -D2,
                                 -0*D2,
                                 mean(dat_for_analysis_4$GH_baseline),
                                 mean(dat_for_analysis_4$higheduce),
                                 mean(dat_for_analysis_4$smoke_baseline)))
  
  
  
  outcome_MSM_3 <- (dat[,]%*%c(1,
                               -D1_plus,
                               -D1_minus,
                               0*D2,
                               D2,
                               mean(dat_for_analysis_4$GH_baseline),
                               mean(dat_for_analysis_4$higheduce),
                               mean(dat_for_analysis_4$smoke_baseline)))
  
  
  
  outcome_MSM_4 <- (dat[,]%*%c(1,
                               -D1_plus,
                               -D1_minus,
                               -0*D2,
                               -D2,
                               mean(dat_for_analysis_4$GH_baseline),
                               mean(dat_for_analysis_4$higheduce),
                               mean(dat_for_analysis_4$smoke_baseline)))
  
  return(list(outcome_MSM_1,
              outcome_MSM_2,
              outcome_MSM_3,
              outcome_MSM_4)) 
}
library(gridExtra)
D_grid <- expand.grid(seq(0,1,0.01),seq(0,1,0.01))
colnames(D_grid) <- c("D1_plus","D1_minus")
#set.seed(36542)
#Z <- rnorm(10000)

dat <- ENGAGE_real_dat_10000_h_3_2_log_10_28_21_alc_coc_0.5_0.25_159[[1]][seq(5000,10000,10),]
var_dat <- ENGAGE_real_dat_10000_h_3_2_log_10_28_21_alc_coc_0.5_0.25_159[[2]][seq(5000,10000,10)]


outcome_list_alc_coc_overall_0 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0,dat,var_dat))

library(doParallel)

#Number of cores/threads
no_cores <- detectCores()-2

cl<-makeCluster(no_cores)
registerDoParallel((cl))
Sys.time()
#clusterExport(cl,varlist=c("ComputeMCBUpperLimitsDesign1"))

upper_credible_limit_list_alc_coc_overall_0 <- parLapply(cl,lapply(outcome_list_alc_coc_overall_0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)
library(latex2exp)



########
outcome_list_alc_coc_overall_0.25 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0.25,dat,var_dat))



#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_overall_0.25 <- parLapply(cl,lapply(outcome_list_alc_coc_overall_0.25,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)







#########
outcome_list_alc_coc_overall_0.5 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0.5,dat,var_dat))


#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_overall_0.5 <- parLapply(cl,lapply(outcome_list_alc_coc_overall_0.5,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)







##############
outcome_list_alc_coc_overall_0.75 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0.75,dat,var_dat))



#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_overall_0.75 <- parLapply(cl,lapply(outcome_list_alc_coc_overall_0.75,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)



#upper_credible_limit_list_alc_coc_overall_0.75<- lapply(lapply(outcome_list_alc_coc_overall_0.75,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)


##########
outcome_list_alc_coc_overall_1.0 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],1.0,dat,var_dat))




#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_overall_1.0 <- parLapply(cl,lapply(outcome_list_alc_coc_overall_1.0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)



#upper_credible_limit_list_alc_coc_overall_1.0<- lapply(lapply(outcome_list_alc_coc_overall_1.0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))




################
################
################
library(gridExtra)
D_grid <- expand.grid(seq(0,1,0.01),seq(0,1,0.01))
colnames(D_grid) <- c("D1_plus","D1_minus")
set.seed(36542)
Z <- rnorm(10000)

dat <- ENGAGE_real_dat_10000_h_2_2_log_10_29_21_TxReadiness_baseline_greater_24_weeks_alc_coc_159[[1]][seq(5000,10000,10),]
var_dat <- ENGAGE_real_dat_10000_h_2_2_log_10_29_21_TxReadiness_baseline_greater_24_weeks_alc_coc_159[[2]][seq(5000,10000,10)]

outcome_list_alc_coc_greater_0 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0,dat,var_dat))

#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_greater_0 <- parLapply(cl,lapply(outcome_list_alc_coc_greater_0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)




#upper_credible_limit_list_alc_coc_greater_0 <- lapply(lapply(outcome_list_alc_coc_greater_0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)



########
outcome_list_alc_coc_greater_0.25 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0.25,dat,var_dat))


#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_greater_0.25 <- parLapply(cl,lapply(outcome_list_alc_coc_greater_0.25,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)



#upper_credible_limit_list_alc_coc_greater_0.25<- lapply(lapply(outcome_list_alc_coc_greater_0.25,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)



#########
outcome_list_alc_coc_greater_0.5 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0.5,dat,var_dat))

#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_greater_0.5 <- parLapply(cl,lapply(outcome_list_alc_coc_greater_0.5,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)




#upper_credible_limit_list_alc_coc_greater_0.5<- lapply(lapply(outcome_list_alc_coc_greater_0.5,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)



##############
outcome_list_alc_coc_greater_0.75 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0.75,dat,var_dat))


#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_greater_0.75 <- parLapply(cl,lapply(outcome_list_alc_coc_greater_0.75,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)



#upper_credible_limit_list_alc_coc_greater_0.75<- lapply(lapply(outcome_list_alc_coc_greater_0.75,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)




##########
outcome_list_alc_coc_greater_1.0 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],1.0,dat,var_dat))


#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_greater_1.0 <- parLapply(cl,lapply(outcome_list_alc_coc_greater_1.0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)



#upper_credible_limit_list_alc_coc_greater_1.0<- lapply(lapply(outcome_list_alc_coc_greater_1.0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)



########
########
########
#######
library(gridExtra)
D_grid <- expand.grid(seq(0,1,0.01),seq(0,1,0.01))
colnames(D_grid) <- c("D1_plus","D1_minus")
set.seed(36542)
Z <- rnorm(10000)

dat <- ENGAGE_real_dat_10000_h_2_2_log_10_28_21_TxReadiness_baseline_lesser_24_weeks_alc_coc_159[[1]][seq(5000,10000,10),]
var_dat <- ENGAGE_real_dat_10000_h_2_2_log_10_28_21_TxReadiness_baseline_lesser_24_weeks_alc_coc_159[[2]][seq(5000,10000,10)]

outcome_list_alc_coc_lesser_0 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0,dat,var_dat))


#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_lesser_0 <- parLapply(cl,lapply(outcome_list_alc_coc_lesser_0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)



#upper_credible_limit_list_alc_coc_lesser_0 <- lapply(lapply(outcome_list_alc_coc_lesser_0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))





########
outcome_list_alc_coc_lesser_0.25 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0.25,dat,var_dat))

#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_lesser_0.25 <- parLapply(cl,lapply(outcome_list_alc_coc_lesser_0.25,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)




#upper_credible_limit_list_alc_coc_lesser_0.25<- lapply(lapply(outcome_list_alc_coc_lesser_0.25,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)




#########
outcome_list_alc_coc_lesser_0.5 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0.5,dat,var_dat))

#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_lesser_0.5 <- parLapply(cl,lapply(outcome_list_alc_coc_lesser_0.5,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)



#upper_credible_limit_list_alc_coc_lesser_0.5<- lapply(lapply(outcome_list_alc_coc_lesser_0.5,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)




##############
outcome_list_alc_coc_lesser_0.75 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0.75,dat,var_dat))


#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_lesser_0.75 <- parLapply(cl,lapply(outcome_list_alc_coc_lesser_0.75,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)



#upper_credible_limit_list_alc_coc_lesser_0.75<- lapply(lapply(outcome_list_alc_coc_lesser_0.75,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))





##########
outcome_list_alc_coc_lesser_1.0 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],1.0,dat,var_dat))


#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_coc_lesser_1.0 <- parLapply(cl,lapply(outcome_list_alc_coc_lesser_1.0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))

library(ggplot2)

stopCluster(cl)


#upper_credible_limit_list_alc_coc_lesser_1.0<- lapply(lapply(outcome_list_alc_coc_lesser_1.0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) colMeans(w))




###########################



g_overall_alc_coc_best_0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_overall_0, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Overall $D_{2,+1} = D_{2,-1} = 0.00$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),
                                                                                                                                                                                                                                                                                                                                                                                                     axis.title.x = element_text(size=22, face="bold"),
                                                                                                                                                                                                                                                                                                                                                                                                     axis.title.y = element_text(size=22, face="bold"),
                                                                                                                                                                                                                                                                                                                                                                                                     axis.text = element_text(size=22,face="bold"))
g_overall_alc_coc_best_0.25 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_overall_0.25, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Overall $D_{2,+1} = D_{2,-1} = 0.25$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))
g_overall_alc_coc_best_0.5 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_overall_0.5, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Overall $D_{2,+1} = D_{2,-1} = 0.50$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))
g_overall_alc_coc_best_0.75 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_overall_0.75, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Overall $D_{2,+1} = D_{2,-1} = 0.75$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))
g_overall_alc_coc_best_1.0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_overall_1.0, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Overall $D_{2,+1} = D_{2,-1} = 1.00$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))


g_greater_alc_coc_best_0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_greater_0, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Greater $D_{2,+1} = D_{2,-1} = 0.00$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))
g_greater_alc_coc_best_0.25 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_greater_0.25, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Greater $D_{2,+1} = D_{2,-1} = 0.25$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))
g_greater_alc_coc_best_0.5 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_greater_0.5, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Greater $D_{2,+1} = D_{2,-1} = 0.50$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))
g_greater_alc_coc_best_0.75 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_greater_0.75, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Greater $D_{2,+1} = D_{2,-1} = 0.75$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))
g_greater_alc_coc_best_1.0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_greater_1.0, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Greater $D_{2,+1} = D_{2,-1} = 1.00$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))

g_lesser_alc_coc_best_0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_lesser_0, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Lesser $D_{2,+1} = D_{2,-1} = 0.00$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))
g_lesser_alc_coc_best_0.25 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_lesser_0.25, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Lesser $D_{2,+1} = D_{2,-1} = 0.25$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))
g_lesser_alc_coc_best_0.5 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_lesser_0.5, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Lesser $D_{2,+1} = D_{2,-1} = 0.50$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))
g_lesser_alc_coc_best_0.75 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_lesser_0.75, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Lesser $D_{2,+1} = D_{2,-1} = 0.75$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'),axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))
g_lesser_alc_coc_best_1.0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'Best EDTR' =as.factor(unlist(lapply(upper_credible_limit_list_alc_coc_lesser_1.0, function(x) (which.min(x)))))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`Best EDTR`)) +geom_tile()+ ggtitle(TeX("Best EDTR Lesser $D_{2,+1} = D_{2,-1} = 1.00$"))+scale_fill_grey()+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(plot.title = element_text(size=22,face="bold"), legend.title = element_text(size=22,face="bold"), legend.text = element_text(size=22,face='bold'), axis.title=element_text(size=22,face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.x = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.title.y = element_text(size=22, face="bold"),                                                                                                                                                                                                                                                                                                                                                                                                      axis.text = element_text(size=22,face="bold"))




pdf("overall_alc_coc_best_0.00.pdf",width=10)
g_overall_alc_coc_best_0
dev.off()

pdf("overall_alc_coc_best_0.25.pdf",width=10)
g_overall_alc_coc_best_0.25
dev.off()

pdf("overall_alc_coc_best_0.50.pdf",width=10)
g_overall_alc_coc_best_0.5

dev.off()

pdf("overall_alc_coc_best_0.75.pdf",width=10)
g_overall_alc_coc_best_0.75
dev.off()

pdf("overall_alc_coc_best_1.00.pdf",width=10)
g_overall_alc_coc_best_1.0
dev.off()


pdf("greater_alc_coc_best_0.00.pdf",width=10)
g_greater_alc_coc_best_0
dev.off()

pdf("greater_alc_coc_best_0.25.pdf",width=10)
g_greater_alc_coc_best_0.25
dev.off()

pdf("greater_alc_coc_best_0.50.pdf",width=10)
g_greater_alc_coc_best_0.5

dev.off()

pdf("greater_alc_coc_best_0.75.pdf",width=10)
g_greater_alc_coc_best_0.75
dev.off()

pdf("greater_alc_coc_best_1.00.pdf",width=10)
g_greater_alc_coc_best_1.0
dev.off()


pdf("lesser_alc_coc_best_0.00.pdf",width=10)
g_lesser_alc_coc_best_0
dev.off()

pdf("lesser_alc_coc_best_0.25.pdf",width=10)
g_lesser_alc_coc_best_0.25
dev.off()

pdf("lesser_alc_coc_best_0.50.pdf",width=10)
g_lesser_alc_coc_best_0.5

dev.off()

pdf("lesser_alc_coc_best_0.75.pdf",width=10)
g_lesser_alc_coc_best_0.75
dev.off()

pdf("lesser_alc_coc_best_1.00.pdf",width=10)
g_lesser_alc_coc_best_1.0
dev.off()

