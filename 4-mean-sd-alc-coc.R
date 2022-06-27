#load("ENGAGE_real_dat_10000_h_3_3_log_7_29_21_12_weeks.rda")


ENGAGE_summary <- cbind(round(colMeans(ENGAGE_real_dat_10000_h_3_2_log_10_28_21_alc_coc_0.5_0.25_159[[1]]),2),
round(apply(ENGAGE_real_dat_10000_h_3_2_log_10_28_21_alc_coc_0.5_0.25_159[[1]],2,sd),2))

round(mean(ENGAGE_real_dat_10000_h_3_2_log_10_28_21_alc_coc_0.5_0.25_159[[2]]),4)

knitr::kable(t(apply(rbind(format(round(colMeans(ENGAGE_real_dat_10000_h_3_2_log_10_28_21_alc_coc_0.5_0.25_159[[1]]),2),nsmall=2),
                round(apply(ENGAGE_real_dat_10000_h_3_2_log_10_28_21_alc_coc_0.5_0.25_159[[1]],2,sd),2)),2,function(x) paste(x[1]," ","(",x[2],")",sep=""))),format="latex")

knitr::kable(rbind(format(round(colMeans(ENGAGE_real_dat_10000_h_3_2_log_10_28_21_alc_coc_0.5_0.25_159[[1]]),2),nsmall=2),
                          round(apply(ENGAGE_real_dat_10000_h_3_2_log_10_28_21_alc_coc_0.5_0.25_159[[1]],2,sd),2)),format="latex")
             
###########

knitr::kable(rbind(format(round(colMeans(ENGAGE_real_dat_10000_h_3_2_log_10_28_21_alc_coc_0.5_0.25_159[[1]]),2),nsmall=2),
                   round(apply(ENGAGE_real_dat_10000_h_3_2_log_10_28_21_alc_coc_0.5_0.25_159[[1]],2,sd),2)),format="latex")

knitr::kable(rbind(format(round(colMeans(ENGAGE_real_dat_10000_h_2_2_log_10_28_21_TxReadiness_baseline_lesser_24_weeks_alc_coc_159[[1]]),2),nsmall=2),
                   round(apply(ENGAGE_real_dat_10000_h_2_2_log_10_28_21_TxReadiness_baseline_lesser_24_weeks_alc_coc_159[[1]],2,sd),2)),format="latex")

knitr::kable(rbind(format(round(colMeans(ENGAGE_real_dat_10000_h_2_2_log_10_29_21_TxReadiness_baseline_greater_24_weeks_alc_coc_159[[1]]),2),nsmall=2),
                   round(apply(ENGAGE_real_dat_10000_h_2_2_log_10_29_21_TxReadiness_baseline_greater_24_weeks_alc_coc_159[[1]],2,sd),2)),format="latex")
################################