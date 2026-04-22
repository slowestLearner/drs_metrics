# --- directly modified from Min's code
rm(list = ls())
library(this.path)
setwd(this.path::this.dir())
library(tidyverse)
library(lubridate)
library(lme4)
library(fixest)
library(dplyr)
library(data.table)
library(tictoc)

# --- perhaps we can match?
load("../../../Data_M/MFdata.RData")
data <- as.data.table(dat)
rm(dat)
data <- data[, .(yyyymm = 100 * year + month, fundid, net_return, gross_return, capm_adj_gret, ff5_adj_gret)]


# match to wficn
tmp <- fread("../../../../Ben-David, Li, Rossi, Song (Investors dont like CAPM)/Data and Code/Morningstar Data/processed/matching/secId_to_wficn.csv")
tmp <- unique(tmp[, .(fundid = FundId, wficn)]) %>% na.omit()
data <- data[tmp, on = .(fundid), nomatch = NULL]
rm(tmp)

# match to my stuff
tmp <- readRDS("../tmp/raw_data/fund_net_ret_risk_adjusted_full_sample.RDS")
tmp <- tmp[, .(yyyymm, wficn, ret, alpha_CAPM, alpha_FF5)]
data <- data[tmp, on = .(yyyymm, wficn), nomatch = NULL]
rm(tmp)

# these match?
data[, cor(net_return, ret)]
data[, cor(capm_adj_gret, alpha_CAPM)]

data[, bin := ntile(capm_adj_gret, 100)]

out <- data[, .(capm_adj_gret = mean(capm_adj_gret), alpha_CAPM = mean(alpha_CAPM)), bin]

ggplot(out, aes(x = capm_adj_gret, y = alpha_CAPM)) +
  geom_point() +
  labs(x = "CAPM-adjusted return (Min)", y = "CAPM-adjusted return (J)") +
  theme(text = element_text(size = 35)) +
  geom_abline(slope = 1, intercept = 0)
