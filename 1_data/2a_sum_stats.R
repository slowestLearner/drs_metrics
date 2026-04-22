# --- Produce summary statistics
library(this.path)
setwd(this.path::this.dir())
source("../runmefirst.R")
options(width = 150)

# Min's data
data <- readRDS("../../../data/funds/min_sample.RDS")

# # numbers to report
# data[, length(unique(fundid))]
# dim(data)[1]

# list of variables
data <- data[order(fundid, yyyymm)]
data[, starting_logtna := first(logtna), fundid]
data[, starting_fund_age := first(fund_age), fundid]

# a bunch of things are in percent
vars_pct <- c(
  "gross_return", "net_return", "benchmark_adj_gret", "capm_adj_gret", "ff3_adj_gret", "carhart_adj_gret", "ff5_adj_gret", "ff6_adj_gret"
  # ,"ner", "turnover"
)
data[, (vars_pct) := lapply(.SD, function(x) x * 100), .SDcols = vars_pct]

vars <- data.table(var = c(
  "gross_return", "net_return", "benchmark_adj_gret",
  "capm_adj_gret",
  "ff3_adj_gret", "carhart_adj_gret", "ff5_adj_gret", "ff6_adj_gret", "tna",
  "lag_logtna", "ner", "fund_age", "turnover", "starting_logtna", "starting_fund_age"
))
vars[, var_idx := .I]

data <- data[, vars[, var], with = F]
data[, idx := .I]
data <- melt(data, id.vars = "idx", variable.name = "var", value.name = "value") %>%
  mutate(var = as.character(var)) %>%
  setDT()
data <- data[vars, on = .(var)]

out <- data[, .(
  obs = .N, obs_non_missing = sum(!is.na(value)),
  x_mean = mean(value, na.rm = T),
  x_sd = sd(value, na.rm = T),
  q01 = quantile(value, 0.01, na.rm = T),
  q05 = quantile(value, 0.05, na.rm = T),
  q25 = quantile(value, 0.25, na.rm = T),
  q50 = quantile(value, 0.50, na.rm = T),
  q75 = quantile(value, 0.75, na.rm = T),
  q95 = quantile(value, 0.95, na.rm = T),
  q99 = quantile(value, 0.99, na.rm = T)
), .(var_idx, var)][order(var_idx)]

to_dir <- "../tmp/raw_data/summary/"
dir.create(to_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(out, paste0(to_dir, "min_sample_sum_stats.RDS"))
