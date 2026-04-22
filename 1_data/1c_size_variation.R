# --- Check the amount of fund size variation. This is only for doing back-of-the-envelope calculations in the paper
library(this.path)
setwd(this.path::this.dir())
source("../utility_functions/runmefirst.R")
options(width = 150)

# data
data <- readRDS("../../data/funds/min_sample.RDS")[, .(yyyymm, fundid, size_1 = exp(lag_logtna))]

# produce a version that scales by monthly fund size changes
tmp <- copy(data)
tmp[, idx := frank(yyyymm, ties.method = "dense")]
tmp <- tmp[tmp[, .(idx = idx - 1, fundid, size = size_1)], on = .(idx, fundid), nomatch = NULL][, idx := NULL]
tmp <- tmp[size != size_1]
tmp[, adj_fund := size / size_1 - 1]
tmp[, adj_fund := Winsorize(adj_fund, quantile(adj_fund, probs = c(.005, .995))), yyyymm]
tmp <- tmp[, .(adj_fund = weighted.mean(adj_fund, size_1)), yyyymm]
tmp <- tmp[order(yyyymm)]

# merge with market return for adjustment
tt <- fread("../../data/factors/ff5_plus_mom_monthly_factors.csv")[, .(yyyymm, mkt = mktRf + rf)]
tmp <- tmp[tt, on = .(yyyymm), nomatch = NULL]
rm(tt)

# get cumulative
tmp <- tmp[order(yyyymm)]
tmp <- tmp[, .(yyyymm,
  cum_adj_fund = cumprod(1 + adj_fund),
  cum_adj_mkt = cumprod(1 + mkt)
)]

data <- data[tmp, on = "yyyymm", nomatch = NULL]
rm(tmp)

data[, size_1_adj_fund := size_1 / cum_adj_fund]
data[, size_1_adj_mkt := size_1 / cum_adj_mkt]
data <- rbindlist(list(
  data[, .(yyyymm, fundid, type = "size", size = size_1)],
  data[, .(yyyymm, fundid, type = "size_adj_fund", size = size_1_adj_fund)],
  data[, .(yyyymm, fundid, type = "size_adj_mkt", size = size_1_adj_mkt)]
))

# Let's figure out fund life length
out <- data[, .(min_size = min(size), max_size = max(size), obs = length(unique(yyyymm))), by = .(fundid, type)]
out[, range := max_size / min_size - 1]
out[, bin := ntile(obs, 20)]

# plot the range of size variation as a function of fund life length
tt <- out[, .(range = median(range), obs = mean(obs / 12)), .(bin, type)]
ggplot(tt[obs <= 20], aes(x = obs, y = range, color = type)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  theme(text = element_text(size = 35)) +
  ggtitle("Size Variation")

# median range of size variation for funds with at least 10 years of data
out[obs >= 10 * 12, median(range), type]
