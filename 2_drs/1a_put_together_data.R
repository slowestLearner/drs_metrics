# --- Put together a regression-ready dataset for DRS estimation
library(this.path)
setwd(this.path::this.dir())
source("../utility_functions/runmefirst.R")
# library(tidyverse)
# library(lubridate)
# library(lme4)
# library(fixest)

# load data
dat <- readRDS("../../../data/funds/min_sample.RDS")[order(fundid, yyyymm)] %>% setDT()
dat$year <- floor(dat$yyyymm / 100)

# get lagged cpi-adjusted tna
dat[, logtna_cpi := log(tna_cpi)]
dat[, idx := frank(yyyymm, ties.method = "dense")]
dat <- merge(dat, dat[, .(idx = idx + 1, fundid, lag_logtna_cpi = log(tna_cpi))], by = c("idx", "fundid"), all.x = T)
dat[is.na(lag_logtna_cpi), lag_logtna_cpi := lag_logtna]

# also get market cap-adjusted
tt <- readRDS("../../../data/macro/total_mkt_cap.RDS")
total_mktcap_base <- tt[yyyymm == max(dat$yyyymm), total_mktcap]
tt[, mktcap_adj := total_mktcap_base / total_mktcap]
tt[, idx := frank(yyyymm, ties.method = "dense")]
tt <- merge(tt[, .(idx, yyyymm, mktcap_adj)], tt[, .(idx = idx + 1, mktcap_adj_1 = mktcap_adj)], by = "idx")[, idx := NULL]
dat <- merge(dat, tt, by = "yyyymm")
rm(tt)

dat[, logtna_mktcap := log(tna * mktcap_adj)]
dat[, lag_logtna_mktcap := log(exp(lag_logtna) * mktcap_adj_1)]
dat[, c("mktcap_adj", "mktcap_adj_1") := NULL]

## fill the missing ner (net expense ratio) with the latest ner information
dat <- dat[order(fundid, yyyymm)]
dat <- dat |>
  group_by(fundid) |>
  fill(ner, .direction = "downup") |>
  ungroup() %>%
  setDT()
dat$ner <- DescTools::Winsorize(dat$ner, quantile(dat$ner, probs = c(0.01, 0.99)))

# require having key variables
dat <- dat[!is.na(expanding_ret)]

# mark if a fund is dead
max_ym <- dat[, max(yyyymm)]
dat <- dat |>
  group_by(fundid) |>
  mutate(dead = max(yyyymm) < max_ym, logage = log(fund_age), si = first(lag_logtna)) %>%
  setDT()
rm(max_ym)

# dat <- dat |> left_join(tmp, by = c("fundid"))
saveRDS(dat, "../tmp/raw_data/reg_table_min.RDS")
