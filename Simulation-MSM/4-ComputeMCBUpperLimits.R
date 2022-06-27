#Compute simultaneous credible intervals from draws of embedded dynamic treatment regime response probabilities given by thetadraws

ComputeMCBUpperLimits <- function(thetadraws,alpha=0.05) {
  
  #Arguments:
  #thetadraws: draws of embedded dynamic treatment regime draws
  #alpha: Probability of excluding the true best EDTR
  
  upper_limit <- rep(NA,4)

  #Compute index of best EDTR
  max_odds_ind <- which.min(colMeans(thetadraws))
  
  #Compute differences between each EDTR and best
  diff_matrix <- matrix(thetadraws[,max_odds_ind],nrow=nrow(thetadraws),ncol=4)-thetadraws

  #Rank differences
  rank_matrix <- apply(diff_matrix[,-max_odds_ind],2,rank,ties.method = 'random')

  #Find max rank
  rank_max <- apply(rank_matrix,1,max)

  #Create sorted differences
  new_dat <- apply(diff_matrix[,],2,sort)

  #Compute 100(1-alpha)% upper quantile
  ranks_quantile <- ceiling(quantile(rank_max,1-alpha))

  #Compute upper limit of credible interval. One for each difference which determines the set of best.
  upper_limit <-new_dat[ranks_quantile,]


  return(upper_limit)
}
