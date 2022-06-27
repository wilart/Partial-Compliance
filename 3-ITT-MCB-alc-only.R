#Create MCB plots showing the optimal and inferior embedded DTRs
library(ggplot2)
library(reshape2)
library(dplyr)

#############
outcome_MSM__1_ITT<-outcome_MSM__2_ITT<-outcome_MSM__3_ITT<-outcome_MSM__4_ITT <- rep(NA,10000)

dat <- dat_for_analysis_4_alc %>% filter(s == 1) %>% rbind(dat_for_analysis_4_alc %>% filter(s == 1) %>% mutate(a2 =-a2), dat_for_analysis_4_alc %>% filter(s==0))

a1 <- dat$a1
s <- dat$s
a2 <- dat$a2
y <- dat$y
GH <- dat$GH_baseline
TR <- dat$higheduce
smoke <- dat$smoke_baseline
weights <- (s==1)/0.5+(s==0)/0.25
beta_matrix <- matrix(NA,nrow=10000,ncol=7)
des_mat <- cbind(1,a1,a2,a1*a2,GH,TR,smoke)
sigma_current <-1
sigma_vec <- rep(1,10000)
set.seed(2364)
for (j in 1:10000) {
  beta_mean <- solve(t(des_mat) %*%diag(weights)%*% des_mat) %*% t(des_mat) %*%diag(weights)%*% y
  beta_var <-  sigma_current * solve(t(des_mat) %*%diag(weights)%*% des_mat) 
  
  
  sigma_current <- LaplacesDemon::rinvchisq(1,159-7,t(y-des_mat%*% beta_mean)%*%(y-des_mat%*% beta_mean)/(159-7))
  sigma_vec[j] <-   sigma_current
  beta_matrix[j,] <- MASS::mvrnorm(1,beta_mean,beta_var)
  print(j)
}

EDTR_ITT <- cbind(apply(beta_matrix,1,function(x) x%*% c(1,1,1,1,mean(GH),mean(TR),mean(smoke))),
                  apply(beta_matrix,1,function(x) x%*% c(1,1,-1,-1,mean(GH),mean(TR),mean(smoke))),
                  apply(beta_matrix,1,function(x) x%*% c(1,-1,1,-1,mean(GH),mean(TR),mean(smoke))),
                  apply(beta_matrix,1,function(x) x%*% c(1,-1,-1,1,mean(GH),mean(TR),mean(smoke))))





ITT_dat <- data.frame(c("EDTR 1",
        "EDTR 2",
        "EDTR 3",
        "EDTR 4"),colMeans(EDTR_ITT))

colnames(ITT_dat) <- c("EDTR","Value")
ggplot(ITT_dat,aes(x=as.factor(EDTR),y=Value,color=as.factor(EDTR))) + geom_point(size=3)

#########
#Compute draws of difference between each embedded DTR and the best
theta_ITT_differences_ <- as.data.frame(t(apply(EDTR_ITT,1,function(x) x[which.min(colMeans(EDTR_ITT))]-x)))
colnames(theta_ITT_differences_) <- c("EDTR 1","EDTR 2", 'EDTR 3',"EDTR 4")

summary_ITT_ <- theta_ITT_differences_%>%melt %>% group_by(variable) %>% summarize(lower=mean(value),Difference=mean(value))
summary_ITT_ <- cbind(summary_ITT_,upper=ComputeMCBUpperLimitsDesign1(EDTR_ITT))


colnames(summary_ITT_) <- c("EDTRs",
                            "lower",
                            "Differences",
                            "upper")
cbPalette= c('#ca0020','#f4a582','#238b45','#404040')
g_plot<- ggplot(summary_ITT_,aes(x=EDTRs,y=Differences,shape=EDTRs,col=EDTRs))+geom_hline(yintercept=0)+geom_pointrange(aes(ymin=lower,ymax=upper),position=position_dodge2(1),size=1.5,fatten=5)+ theme_classic()+theme(
  axis.title.x = element_text(size=14, face="bold"),
  axis.title.y = element_text(size=14, face="bold"),
  axis.text = element_text(size=14,face="bold"),
  legend.position="none"
)+  scale_shape_manual(values = (c("EDTR 1" = 15,"EDTR 2" = 16,"EDTR 3" = 17,"EDTR 4" = 18)) )+scale_colour_manual(values =(c("EDTR 1" = cbPalette[1],"EDTR 2" =cbPalette[2],"EDTR 3" = cbPalette[3],"EDTR 4" = cbPalette[4]))  )+ggtitle("MCB Plot ENGAGE Overall Alcohol Only")+scale_y_continuous(limits = c(-6.5,5))+theme(plot.title = element_text(size=22))
g_plot

#pdf("Real-data-MCB-plot-overall-24-weeks-3-10-29-21-alc-only.pdf",width=8,height=8)
g_plot
#dev.off()
##############################
########################
##############################