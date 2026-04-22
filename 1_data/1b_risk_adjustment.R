# --- Compute factor-adjusted fund returns. Do both rolling and full-sample.

# load libraries
library(this.path)
setwd(this.path::this.dir())
source("../utility_functions/runmefirst.R")

# get each fund, use rolling 36m regressions to compute alphas
data <- readRDS("../../data/funds/active_equity_fund_monthly_data.RDS")[, .(yyyymm, wficn, ret)]
data[, ret := Winsorize(ret, quantile(ret, probs = c(.005, .995))), yyyymm]

# merge with factors and take out rf
tmp <- fread("../../data/factors/ff5_plus_mom_monthly_factors.csv")
data <- data[tmp, on = .(yyyymm), nomatch = NULL]
rm(tmp)
data[, ret := (1 + ret) / (1 + rf) - 1]
data[, rf := NULL]

# model specifications
spec_data <- data.table(
  model_name = c("CAPM", "FF3", "FFC", "FF5", "FF6"),
  formula = c(
    "ret ~ mktRf",
    "ret ~ mktRf + smb + hml",
    "ret ~ mktRf + smb + hml + mom",
    "ret ~ mktRf + smb + hml + rmw + cma",
    "ret ~ mktRf + smb + hml + rmw + cma + mom"
  )
)
spec_data[, spec_idx := .I]

data[, idx := frank(yyyymm, ties.method = "dense")]
data_all <- copy(data)
rm(data)

# ---- 1) estimate using rolling window

lookback_period <- 36
min_obs <- 24
time_data <- unique(data_all[, .(idx, yyyymm)])[order(idx)]
time_data <- time_data[idx >= (lookback_period + 1)]

to_dir <- "../tmp/raw_data/fund_net_ret_risk_adjusted/"
dir.create(to_dir, recursive = T, showWarnings = F)

# take out those already done
time_data <- time_data[!(yyyymm %in% as.integer(gsub(".RDS", "", list.files(to_dir))))]

for (this_idx in time_data[, idx]) {
  tic(time_data[idx == this_idx, yyyymm])
  # 1. Filter the rolling window for the whole month
  data_month <- data_all[idx %in% (this_idx - lookback_period):this_idx]

  # 2. Perform the regressions within the data.table using 'by'
  out <- data_month[,
    {
      # Check if this specific fund has enough observations in the train period
      train_data <- .SD[idx < this_idx]
      test_data <- .SD[idx == this_idx]

      if (nrow(train_data) >= min_obs && nrow(test_data) > 0) {
        # List to store alphas for this fund
        results <- list(ret = test_data$ret)

        # Loop through models defined in your spec_data
        for (i in seq_len(nrow(spec_data))) {
          m_name <- spec_data$model_name[i]
          m_form <- as.formula(spec_data$formula[i])

          # Fit model on training data
          mm <- lm(m_form, data = train_data)

          # Calculate alpha: actual ret - predicted ret
          # Use [[1]] because test_data is a single month
          alpha_val <- test_data$ret - predict(mm, newdata = test_data)
          results[[paste0("alpha_", m_name)]] <- as.numeric(alpha_val)
        }
        results # This is returned for each wficn
      }
    },
    by = wficn
  ] # This replaces the split() and future_lapply()

  # Save the month's results
  saveRDS(out, paste0(to_dir, time_data[idx == this_idx, yyyymm], ".RDS"))
  toc()
}

# combine together. some of them forgot to include yyyymm
p.get_one <- function(this_file) {
  this_ym <- as.integer(substr(this_file, 1, 6))
  data <- readRDS(paste0(to_dir, this_file))
  data[, yyyymm := this_ym]
  return(data)
}

files <- list.files(to_dir)

tic()
out <- rbindlist(lapply(files, p.get_one), use.names = T)
toc()

# combine into a single file
to_file <- paste0(substr(to_dir, 1, nchar(to_dir) - 1), ".RDS")
saveRDS(out, to_file)
unlist(to_dir, recursive = T)

# ---- 2) estimate using full sample. Require at least 24 obs


# require some min num of obs
min_obs <- 24
data <- copy(data_all)
data[, obs := .N, wficn]
data <- data[obs >= min_obs][, obs := NULL]


# Perform the regressions within the data.table using 'by'
tic()
out <- data[,
  {
    results <- list(yyyymm = .SD$yyyymm, ret = .SD$ret)
    # Loop through models defined in your spec_data
    for (i in seq_len(nrow(spec_data))) {
      m_name <- spec_data$model_name[i]
      m_form <- as.formula(spec_data$formula[i])

      # Fit model on training data
      mm <- lm(m_form, data = .SD)

      # Calculate alpha: actual ret - predicted ret
      alpha_val <- .SD$ret - predict(mm)
      results[[paste0("alpha_", m_name)]] <- as.numeric(alpha_val)
    }
    results
  },
  by = wficn
]
toc()

saveRDS(out, "../tmp/raw_data/fund_net_ret_risk_adjusted_full_sample.RDS")
