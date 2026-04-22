# DRS Metrics Code Documentation



## `1_data/`
basic data processing

- `1_data/1a_move_data.R`: keeps track of where some key input data comes from. 
- `1_data/1a1_min_data.R`: Loads legacy `MFdata.RData` and saves as data.table
- `1_data/1b_risk_adjustment.R`: Computes fund-level factor-adjusted returns using both rolling 36-month and full-sample regressions.
- `1_data/1c_size_variation.R`: Quantifies within-fund size variation over fund lifetimes. Just for back-of-the-envelope reporting in the paper.
- `1_data/2a_sum_stats.R`: Produces summary statistics and save in a table. 

source("~/.runmefirst")

## `2_drs/`

- `2_drs/1a_put_together_data.R`: Prepares a regression-ready dataset. 
- `2_drs/1b_existing_methods.R`: Runs only the "existing methods" DRS estimators (OLS/FE/RD/RE). 
- `2_drs/1c_wb.R`: Main WB estimation for average DRS. 
- `2_drs/1d_wb_heterogeneity.R`: Estimates WB heterogeneity by fund-level tracking error, expense ratio, and SMB loading bins across multiple return definitions.
- `2_drs/1d_wb_heterogeneity_trying.R`: Exploratory baseline heterogeneity checks using alternative binning choices and covariates (turnover, tracking error, fees, SMB, liquidity, size).
- 

## `3a_simulation/`

- `3a_simulation/1_min_script_archive.R`: Archived Monte Carlo generator and estimator benchmarking FE/RD/RE/WB/CW under varying DRS and selection scenarios.
- `3a_simulation/1a_calibrate_survival.R`: Calibrates empirical cumulative fund-death probabilities by horizon (1-20 years) from the Min panel for simulation inputs.
- `3a_simulation/1b_min_script.R`: Simulates fund panels with correlated skill/size dynamics and compares estimator bias, SD, RMSE, and rejection rates across exit rules.
- `3a_simulation/1c_min_script_longer_history.R`: Extends the simulation design to longer histories and custom survival rates to test estimator behavior in long panels.
- `3a_simulation/1d_updated_wb.R`: Tests a modified "new WB" estimator against FE/RD/RE/WB in the simulation environment and exports comparative performance metrics.
- `3a_simulation/2a_check_results.R`: Aggregates simulation output files and plots how RE bias changes across alternative parameterizations.
- `3a_simulation/2b_re_zero_point.R`: Compares simulation variants to diagnose conditions where RE bias approaches zero (notably under altered size-law coefficients).
- `3a_simulation/2c_check_longer_sample.R`: Visualizes estimator bias versus `cor(a_i, phi_i)` for 240-period simulations without selected size-process channels.
- `3a_simulation/2d_new_method.R`: Plots and compares FE/RD/RE/WB/new-WB bias profiles from updated 240-period simulation outputs by exit regime.

## `3b_simulation/`

- `3b_simulation/1_drs.R`: Parallelized large-scale DRS simulation runner that sweeps structural parameters and writes per-specification result objects.
- `3b_simulation/2_drs_summary.R`: Reads simulation outputs and generates formatted summary/figure-ready comparisons across methods and fund-exit mechanisms.
- `3b_simulation/2_drs_summary_achive.R`: Archived plotting/summary script for earlier simulation output layouts and figure templates.

## `4_one_off/`

- `4_one_off/1a_alpha_inputs.R`: Constructs portfolio expected returns, covariance matrices, and benchmark weights (full and split samples) for constrained active optimization.
- `4_one_off/1b_alpha_optimization.R`: Solves short-constrained active portfolio QPs with OSQP, calibrates risk aversion to target active vol, and evaluates OOS active performance.
- `4_one_off/1b_alpha_optimization_archive.R`: Archived in-sample-only version of the OSQP active portfolio optimization workflow.
- `4_one_off/1c_alpha_visualization_todel.R`: Loads optimization outputs, annualizes IS/OOS active moments, and merges portfolio metadata for visualization diagnostics.
- `4_one_off/2a_rd_quantification.R`: Uses tabulated Zhu (2018) size-effect estimates to plot linear/log specifications and numerically integrate total decline.
- `4_one_off/3_mechanism_illustration.R`: Creates a toy time-series illustration comparing raw versus forward-demeaned returns to explain the RD-style transformation mechanism.
- `4_one_off/4_performance_persistence.R`: Estimates and visualizes yearly size-growth sensitivity to past gross/benchmark-adjusted fund performance.

## `produce_figures/`

- `produce_figures/simulation_stylized.R`: Produces stylized simulation plots (bias/SD/RMSE versus `cor(a_i, phi_i)`) by estimator and exit type for presentation figures.
- `produce_figures/wb_estimates.R`: Visualizes WB coefficient estimates with confidence intervals across return definitions and size-scaling specifications.

## `utility_functions/`

- `utility_functions/runmefirst.R`: load libraries
- `utility_functions/Analysis_DRS_forJ.R`: Comprehensive legacy analysis script defining core DRS/WB estimators and running baseline plus heterogeneity empirical analyses.
- `utility_functions/cluster2.R`: Implements one-way/two-way clustered covariance and test utilities (`coeftest.cluster`, `coeftest.cluster.re`) used across estimation scripts.
- `utility_functions/cluster2_archive.R`: Archived version of the clustered standard-error helper functions with older formatting/logic.
- `utility_functions/updated_functions.R`: Revised robust-RE clustering helpers that align model-used rows explicitly to avoid misalignment in random-effect residual adjustments.