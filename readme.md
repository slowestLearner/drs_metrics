# DRS Metrics Code Documentation



## `1_data/`
basic data processing

- `1_data/1a_move_data.R`: keeps track of where some key input data comes from. 
- `1_data/1a1_min_data.R`: Loads legacy `MFdata.RData` and saves as data.table
- `1_data/1b_risk_adjustment.R`: Computes fund-level factor-adjusted returns using both rolling 36-month and full-sample regressions.
- `1_data/1c_size_variation.R`: Quantifies within-fund size variation over fund lifetimes. Just for back-of-the-envelope reporting in the paper.
- `1_data/2a_sum_stats.R`: Produces summary statistics and save in a table. 


## `2_drs/`

- `2_drs/1a_put_together_data.R`: Prepares a regression-ready dataset. 
- `2_drs/1b_existing_methods.R`: Runs only the "existing methods" DRS estimators (OLS/FE/RD/RE). 
- `2_drs/1c_wb.R`: Main WB estimation for average DRS. 
- `2_drs/1d_wb_heterogeneity.R`: Estimates WB heterogeneity by fund-level tracking error, expense ratio, and SMB loading bins across multiple return definitions.
- `2_drs/1d_wb_heterogeneity_trying.R`: Exploratory baseline heterogeneity checks using alternative binning choices and covariates (turnover, tracking error, fees, SMB, liquidity, size). 



## `3_simulation/`

- `3_simulation/1a_calibrate_survival.R`: Calibrates empirical cumulative fund survival probabilities by horizon (1-20 years). 
- `3_simulation/1b_wb.R`: Tests a modified "new WB" estimator against FE/RD/RE/WB in the simulation environment and exports comparative performance metrics.

## `4_one_off/`

- `4_one_off/1a_alpha_inputs.R`: Constructs portfolio expected returns, covariance matrices, and benchmark weights (full and split samples) for constrained active optimization.
- `4_one_off/1b_alpha_optimization.R`: Solves short-constrained active portfolio QPs with OSQP, calibrates risk aversion to target active vol, and evaluates OOS active performance.
- `4_one_off/2a_rd_quantification.R`: Uses tabulated Zhu (2018) size-effect estimates to plot linear/log specifications and numerically integrate total decline.
- `4_one_off/3_mechanism_illustration.R`: Creates a toy time-series illustration comparing raw versus forward-demeaned returns to explain the RD-style transformation mechanism.
- `4_one_off/4_performance_persistence.R`: Estimates yearly size-growth sensitivity. Used in the paper for a back-of-the-envelope calculation

## `produce_figures/`

- `produce_figures/simulation_stylized.R`: Produces stylized simulation plots (bias/SD/RMSE versus `cor(a_i, phi_i)`) by estimator and exit type for presentation figures.
- `produce_figures/wb_estimates.R`: Visualizes WB coefficient estimates with confidence intervals across return definitions and size-scaling specifications.

## `utility_functions/`

- `utility_functions/runmefirst.R`: load libraries
- `utility_functions/cluster2_archive.R`: Implements one-way/two-way clustered covariance and test utilities (`coeftest.cluster`, `coeftest.cluster.re`) used across estimation scripts.
