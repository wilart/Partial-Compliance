Simulation:
"1-Sim-Full-MSM-beta.R": Scenario 2 simulation
"1-Sim-Full-MSM-truncnorm-deciles.R" Scenario 1 simulaton with discrete compliance
"1-Sim-Full-MSM-truncnorm-increase-response-rate.R" Scenario 1 with larger engagement/response rate.
"1-Sim-Full-MSM-truncnorm-zero-inflated-160.R" Zero-inflated simulation
"1-Sim-Full-MSM-truncnorm.R" Scenario 1 simulation


Run method
"MSM-run-Scenario-1-TN-200.R" Run scenario 1 with 200 sample size
"MSM-run-Scenario-1-TN-400.R" Run scenario 1 with 400 sample size
"MSM-run-Scenario-1-TN-deciles.R" Run Scenario 1 with discrete compliance
"MSM-run-Scenario-1-TN-higher-response-rate-200.R" Run Scenario 1 with higher engagement/response rate
"MSM-run-Scenario-2-Beta-no-covariates-200.R" Run Scenario 2 with 200 sample size
"MSM-run-Scenario-2-Beta-no-covariates-400.R" Run Scenario 2 with 400 sample size
"MSM-run-Zero-Inflation-TN-160.R"  Run Zero-inflated simulation mimicking data

Fit model
"1-dtmvnorm.R" Evaluating multivariate truncated normal density
"2-Beta-fit-no-covariates.R" Fit model with covariates in outcome model but not compliance model
"2-Beta-fit-with-covariates-2.R" Fit model with covariates in outcome and compliance model

Multiple comparisons with the best
"4-ComputeMCBUpperLimits.R" Compute upper credible interval limits using Bayesian multiple comparisons with the best
"5-Compliance-coverage-computation.R" Compute Bayesian coverage of credible intervals based off multiple comparisons with the best and ignoring multiplicity

Calculate Bias/SD
"3-Bias-SD-Scenario-1-TN.R" Bias/SD for Scenario 1
"3-Bias-SD-Scenario-2-Beta.R" Bias/SD for Scenario 2
"3-Bias-SD-TN-higher-response-rate.R" Bias/SD for Scenario 1 with increased engagement/response rate
"3-Bias-SD-with-Scenario-1-TN-deciles.R" Bias/SD for Scenario 1 with discrete potential compliance
"3-Bias-SD-Zero-Inflation-TN-160.R" Bias/SD for Zero Inflation with sample size of 160

Real data analysis:
"Synthetic-data" Simulates data similar to the ENGAGE study.
"1-dtmvnorm-for-TN-Scenario-1.R" Compute density of multivariate truncated normal
"2-Beta-fit-with-covariates-2.R" Fit with covariates in both outcome and compliance model
"2-real-data-MSM-run-alc-coc.R" Run Gibbs sampler for alcohol + cocaine outcome
"2-real-data-MSM-run-alc-only.R" Run Gibbs sampler for alcohol only outcome
"3-ITT-MCB-alc-coc.R" Make plot for upper credible intervals for intention-to-treat analysis for alcohol + cocaine outcome
"3-ITT-MCB-alc-only.R" Make plot for upper credible intervals for intention-to-treat analysis for alcohol only outcome
"4-mean-sd-alc-coc.R" Compute mean and SD for alcohol+cocaine outcome
"5-heatmap-alc-coc-by-D2-best.R" Create heatmaps for alcohol+cocaine outcome set of best membership
"5-heatmap-alc-only-by-D2-best.R" Create heatmaps for alcohol only outcome set of best membership
"6-alc_coc_histogram_outcome.R" Create histogram of transformed alcohol + cocaine outcome.



