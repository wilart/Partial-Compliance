

set.seed(3275)
ENGAGE_real_dat_10000_h_3_2_log_10_28_21_alc_coc_0.5_0.25_159 <- computeBetas(10000,dat_for_analysis_4,H_responders = 3,H_non_responders = 2)



#subgroup analysis


dat_for_analysis_subgroup_greater <- dat_for_analysis_4 %>% filter(TxReadiness_baseline>=median(TxReadiness_baseline))


set.seed(3758512)
ENGAGE_real_dat_10000_h_2_2_log_10_29_21_TxReadiness_baseline_greater_24_weeks_alc_coc_159 <- computeBetas(10000,dat_for_analysis_subgroup_greater,H_responders = 2,H_non_responders = 2)

dat_for_analysis_subgroup_less <- dat_for_analysis_4 %>% filter(TxReadiness_baseline<median(TxReadiness_baseline))

set.seed(3758561)
ENGAGE_real_dat_10000_h_2_2_log_10_28_21_TxReadiness_baseline_lesser_24_weeks_alc_coc_159 <- computeBetas(10000,dat_for_analysis_subgroup_less,H_responders = 2,H_non_responders = 2)
