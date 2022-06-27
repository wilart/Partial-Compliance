complianceMatrix <- function(D1_plus,D1_minus,D2,dat,var_dat) {
  
  outcome_MSM_1 <- outcome_MSM_2 <- outcome_MSM_3 <- outcome_MSM_4 <- rep(NA,10000)
  
  outcome_MSM_1 <- (dat[,]%*%c(1,
                               D1_plus,
                               D1_minus,
                               D2,
                               0*D2,
                               mean(dat_for_analysis_4_alc$GH_baseline),
                               mean(dat_for_analysis_4_alc$higheduce),
                               mean(dat_for_analysis_4_alc$smoke_baseline)))
  
  
  
  outcome_MSM_2 <- (dat[,]%*%  c(1,
                                 D1_plus,
                                 D1_minus,
                                 -D2,
                                 -0*D2,
                                 mean(dat_for_analysis_4_alc$GH_baseline),
                                 mean(dat_for_analysis_4_alc$higheduce),
                                 mean(dat_for_analysis_4_alc$smoke_baseline)))
  
  
  
  outcome_MSM_3 <- (dat[,]%*%c(1,
                               -D1_plus,
                               -D1_minus,
                               0*D2,
                               D2,
                               mean(dat_for_analysis_4_alc$GH_baseline),
                               mean(dat_for_analysis_4_alc$higheduce),
                               mean(dat_for_analysis_4_alc$smoke_baseline)))
  
  
  
  outcome_MSM_4 <- (dat[,]%*%c(1,
                               -D1_plus,
                               -D1_minus,
                               -0*D2,
                               -D2,
                               mean(dat_for_analysis_4_alc$GH_baseline),
                               mean(dat_for_analysis_4_alc$higheduce),
                               mean(dat_for_analysis_4_alc$smoke_baseline)))
  
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

dat <- ENGAGE_real_dat_10000_h_3_2_log_10_29_21_alc_only_0.5_0.25_159[[1]][seq(5000,10000,10),]
var_dat <- ENGAGE_real_dat_10000_h_3_2_log_10_29_21_alc_only_0.5_0.25_159[[2]][seq(5000,10000,10)]


outcome_list_alc_only_overall_0 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0,dat,var_dat))

library(doParallel)

#Number of cores/threads
no_cores <- detectCores()-2

cl<-makeCluster(no_cores)
registerDoParallel((cl))
Sys.time()
clusterExport(cl,varlist=c("ComputeMCBUpperLimitsDesign1"))

upper_credible_limit_list_alc_only_overall_0 <- parLapply(cl,lapply(outcome_list_alc_only_overall_0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) ComputeMCBUpperLimitsDesign1(w))

library(ggplot2)
library(latex2exp)



########
outcome_list_alc_only_overall_0.25 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0.25,dat,var_dat))



#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_only_overall_0.25 <- parLapply(cl,lapply(outcome_list_alc_only_overall_0.25,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) ComputeMCBUpperLimitsDesign1(w))

library(ggplot2)







#########
outcome_list_alc_only_overall_0.5 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0.5,dat,var_dat))


#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_only_overall_0.5 <- parLapply(cl,lapply(outcome_list_alc_only_overall_0.5,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) ComputeMCBUpperLimitsDesign1(w))

library(ggplot2)







##############
outcome_list_alc_only_overall_0.75 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],0.75,dat,var_dat))



#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_only_overall_0.75 <- parLapply(cl,lapply(outcome_list_alc_only_overall_0.75,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) ComputeMCBUpperLimitsDesign1(w))

library(ggplot2)



#upper_credible_limit_list_alc_only_overall_0.75<- lapply(lapply(outcome_list_alc_only_overall_0.75,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) ComputeMCBUpperLimitsDesign1(w))

library(ggplot2)


##########
outcome_list_alc_only_overall_1.0 <-apply(D_grid,1,function(x) complianceMatrix(x[1],x[2],1.0,dat,var_dat))




#Number of cores/threads




Sys.time()


upper_credible_limit_list_alc_only_overall_1.0 <- parLapply(cl,lapply(outcome_list_alc_only_overall_1.0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) ComputeMCBUpperLimitsDesign1(w))

library(ggplot2)



#upper_credible_limit_list_alc_only_overall_1.0<- lapply(lapply(outcome_list_alc_only_overall_1.0,function(x) cbind(x[[1]][],x[[2]][],x[[3]][],x[[4]][])),function(w) ComputeMCBUpperLimitsDesign1(w))




################
################
################
library(gridExtra)


###########################
g1_overall_alc_only_0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 1\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0, function(x) x[1]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 1\n Included`)) +geom_tile()+ ggtitle(TeX("1st EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g2_overall_alc_only_0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 2\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0, function(x) x[2]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 2\n Included`)) +geom_tile()+ ggtitle(TeX("2nd EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g3_overall_alc_only_0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 3\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0, function(x) x[3]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 3\n Included`)) +geom_tile()+ ggtitle(TeX("3rd EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g4_overall_alc_only_0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 4\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0, function(x) x[4]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 4\n Included`)) +geom_tile() + ggtitle(TeX("4th EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))

g1_overall_alc_only_0.25 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 1\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0.25, function(x) x[1]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 1\n Included`)) +geom_tile()+ ggtitle(TeX("1st EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g2_overall_alc_only_0.25 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 2\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0.25, function(x) x[2]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 2\n Included`)) +geom_tile()+ ggtitle(TeX("2nd EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g3_overall_alc_only_0.25 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 3\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0.25, function(x) x[3]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 3\n Included`)) +geom_tile()+ ggtitle(TeX("3rd EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g4_overall_alc_only_0.25 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 4\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0.25, function(x) x[4]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 4\n Included`)) +geom_tile() + ggtitle(TeX("4th EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))

g1_overall_alc_only_0.5 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 1\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0.5, function(x) x[1]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 1\n Included`)) +geom_tile()+ ggtitle(TeX("1st EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g2_overall_alc_only_0.5 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 2\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0.5, function(x) x[2]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 2\n Included`)) +geom_tile()+ ggtitle(TeX("2nd EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g3_overall_alc_only_0.5 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 3\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0.5, function(x) x[3]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 3\n Included`)) +geom_tile()+ ggtitle(TeX("3rd EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g4_overall_alc_only_0.5 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 4\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0.5, function(x) x[4]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 4\n Included`)) +geom_tile() + ggtitle(TeX("4th EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))

g1_overall_alc_only_0.75 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 1\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0.75, function(x) x[1]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 1\n Included`)) +geom_tile()+ ggtitle(TeX("1st EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g2_overall_alc_only_0.75 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 2\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0.75, function(x) x[2]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 2\n Included`)) +geom_tile()+ ggtitle(TeX("2nd EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g3_overall_alc_only_0.75 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 3\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0.75, function(x) x[3]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 3\n Included`)) +geom_tile()+ ggtitle(TeX("3rd EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g4_overall_alc_only_0.75 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 4\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_0.75, function(x) x[4]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 4\n Included`)) +geom_tile() + ggtitle(TeX("4th EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))

g1_overall_alc_only_1.0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 1\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_1.0, function(x) x[1]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 1\n Included`)) +geom_tile()+ ggtitle(TeX("1st EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g2_overall_alc_only_1.0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 2\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_1.0, function(x) x[2]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 2\n Included`)) +geom_tile()+ ggtitle(TeX("2nd EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g3_overall_alc_only_1.0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 3\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_1.0, function(x) x[3]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 3\n Included`)) +geom_tile()+ ggtitle(TeX("3rd EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))
g4_overall_alc_only_1.0 <- ggplot2::ggplot(data.frame(cbind(D_grid,'EDTR 4\n Included' =unlist(lapply(upper_credible_limit_list_alc_only_overall_1.0, function(x) x[4]>=0))),check.names = F),aes(x=D1_plus,y=D1_minus,fill=`EDTR 4\n Included`)) +geom_tile() + ggtitle(TeX("4th EDTR Alcohol Overall"))+scale_fill_manual(values=c("FALSE" = "grey","TRUE"="black"))+xlab(TeX("D_{1,+1}")) + ylab(TeX("D_{1,-1}"))+theme(axis.title=element_text(size=14,face="bold"))


library(gridExtra)
library(grid)


pdf("overall_alc_only_0.00.pdf",width=10)
grid.arrange(g1_overall_alc_only_0,
             g2_overall_alc_only_0,
             g3_overall_alc_only_0,
             g4_overall_alc_only_0,nrow=2,top=textGrob(TeX("Inclusion in Set of Best $D_{2,+1} = D_{2,-1} = 0.00$"),gp = gpar(fontsize = 20, col = 'black', fontface = 'bold')))
dev.off()

pdf("overall_alc_only_0.25.pdf",width=10)
grid.arrange(g1_overall_alc_only_0.25,
             g2_overall_alc_only_0.25,
             g3_overall_alc_only_0.25,
             g4_overall_alc_only_0.25,nrow=2,top=textGrob(TeX("Inclusion in Set of Best $D_{2,+1} = D_{2,-1} = 0.25$"),gp = gpar(fontsize = 20, col = 'black', fontface = 'bold')))
dev.off()

pdf("overall_alc_only_0.50.pdf",width=10)

grid.arrange(g1_overall_alc_only_0.5,
             g2_overall_alc_only_0.5,
             g3_overall_alc_only_0.5,
             g4_overall_alc_only_0.5,nrow=2,top=textGrob(TeX("Inclusion in Set of Best $D_{2,+1} = D_{2,-1} = 0.50$"),gp = gpar(fontsize = 20, col = 'black', fontface = 'bold')))
dev.off()

pdf("overall_alc_only_0.75.pdf",width=10)

grid.arrange(g1_overall_alc_only_0.75,
             g2_overall_alc_only_0.75,
             g3_overall_alc_only_0.75,
             g4_overall_alc_only_0.75,nrow=2,top=textGrob(TeX("Inclusion in Set of Best $D_{2,+1} = D_{2,-1} = 0.75$"),gp = gpar(fontsize = 20, col = 'black', fontface = 'bold')))
dev.off()

pdf("overall_alc_only_1.00.pdf",width=10)

grid.arrange(g1_overall_alc_only_1.0,
             g2_overall_alc_only_1.0,
             g3_overall_alc_only_1.0,
             g4_overall_alc_only_1.0,nrow=2,top=textGrob(TeX("Inclusion in Set of Best $D_{2,+1} = D_{2,-1} = 1.00$"),gp = gpar(fontsize = 20, col = 'black', fontface = 'bold')))
dev.off()

