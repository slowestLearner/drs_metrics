# --- i think if possible, we should get portfolio-weighted big-ask spreads
# when bid-ask is missing, use market cap-imputed
# if still missing, just use SMB loadings
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



# load some base functions (for clustering)
source("../../../data_M/cluster2.R")

# function to compute robust standard errors
robust.se.remodel <- function(data, model, cluster = c("wficn", "yyyymm"), weights = FALSE) {
  dd <- data.frame(y = model.frame(model)[, 1], model.matrix(model)[, -1])

  random_intercepts <- ranef(model)$wficn
  dd$raneff <- rep(random_intercepts[, 1], as.numeric(table(data$wficn)))

  # subtract random intercepts
  dd$y <- dd$y - dd$raneff
  xnames <- setdiff(names(dd)[-1], "raneff")
  fmla <- paste0("y ~ ", paste(xnames, collapse = " + "))
  if (weights) {
    m1 <- lm(fmla, data = dd, weights = w)
  } else {
    m1 <- lm(fmla, data = dd)
  }

  dd$wficn <- data$wficn
  dd$yyyymm <- data$yyyymm
  dd$yyyy <- data$yyyy
  dd$year <- data$year

  # model.frame(model)  # Get the actual data used in fitting with the 1st column always the response variable
  # model.matrix(model)   # Get the model matrix (with interaction terms expanded)

  if (length(cluster) == 1) {
    return(coeftest.cluster.re(data = dd, fm = m1, cluster1 = cluster, weights = weights))
  }
  if (length(cluster) == 2) {
    return(coeftest.cluster.re(data = dd, fm = m1, cluster1 = cluster[1], cluster2 = cluster[2], weights = weights))
    # if cluster is c("wficn", "yyyymm"), then we need to cluster by both wficn and yyyymm
  }
}


#######################################################
#### Load data
#######################################################

# data: merge with fund size, tc
data <- readRDS("../tmp/raw_data/fund_net_ret_risk_adjusted.RDS")
tmp <- readRDS("../../../data/funds/active_equity_fund_monthly_data.RDS")[, .(yyyymm, wficn, exp_ratio, lag_logtna = log(mtna_1), turnover = turn_ratio)]
tmp[, idx := frank(yyyymm, ties.method = "dense")]
tmp <- tmp[tmp[, .(idx = idx - 1, wficn, logtna = lag_logtna)], on = .(idx, wficn), nomatch = NULL][, idx := NULL]
data <- data[tmp, on = .(yyyymm, wficn), nomatch = NULL]
data <- data[exp_ratio < .05] # otherwise a bit weird
rm(tmp)

# sometimes cluster by year instead of year-month (seems to make little difference)
data[, yyyy := floor(yyyymm / 100)]
data[, year := floor(yyyymm / 100)]

# infer fund age from sample start
tmp <- unique(data[, .(yyyymm)])
tmp[, date := as.Date(as.character(yyyymm * 100 + 28), "%Y%m%d")]
data <- data[tmp, on = .(yyyymm)]
data[, min_date := min(date), by = wficn]
data[, fund_age := as.numeric(date - min_date) / 365][, min_date := NULL]
data[, logage := log(fund_age + 1)]
rm(tmp)

# add exp_ratio back into performance measures. Also turn these variables into percent
vars <- names(data)
vars <- vars[grepl("alpha|ret", vars)]
for (this_var in vars) {
  data[, (this_var) := (get(this_var) + exp_ratio / 12) * 100]
}
rm(vars, this_var)

# sort
data <- data[order(yyyymm, wficn)]

# mark if a fund is dead by the end of the sample
max_yyyymm <- max(data$yyyymm)
data <- data |>
  group_by(wficn) |>
  mutate(dead = max(yyyymm) < max_yyyymm, logage = log(fund_age + 1), si = first(lag_logtna)) %>%
  setDT()
rm(max_yyyymm)
data_bk <- copy(data)

### estimate size process (different fixed effect choices)

sizefe1 <- feols(logtna ~ lag_logtna + ret + logage + exp_ratio | wficn + yyyymm, cluster = c("wficn", "yyyymm"), data = data)
sizefe2 <- feols(logtna ~ lag_logtna + ret + logage + exp_ratio | wficn, cluster = c("wficn", "yyyymm"), data = data)

## Table 2: size dynamics
summary(sizefe1)
summary(sizefe2)

## within-between estimation
wb <- function(femod, data = data, ret.name = "alpha_CAPM", time_cluster = "yyyymm") {
  int <- fixef(femod)$wficn
  data$phi <- rep(int, as.numeric(table(data$wficn)))
  mod <- lmer(as.formula(paste0(ret.name, " ~ lag_logtna + phi + si + (1 | wficn)")), data = data)
  obj <- robust.se.remodel(data = data, model = mod, cluster = c("wficn", time_cluster))
  return(obj)
}

# different size processes
wb(sizefe1, data, "alpha_CAPM")
wb(sizefe1, data, "alpha_FF3")


######################### fullsample, two-way interaction of size and portfolio characteristics

# # demonstration: these things change over time
# out <- data[, .(
#   turnover = mean(turnover), exp_ratio = mean(exp_ratio),
#   r2_ff6 = var(alpha_FF6) / var(ret)
# ), year]

# ggplot(out, aes(x = year, y = turnover)) +
#   geom_line(lwd = 1) +
#   geom_hline(yintercept = 0, lty = 3) +
#   ggtitle("turnover") +
#   theme_classic() +
#   theme(text = element_text(size = 35))

# ggplot(out, aes(x = year, y = exp_ratio)) +
#   geom_line(lwd = 1) +
#   geom_hline(yintercept = 0, lty = 3) +
#   ggtitle("exp_ratio") +
#   theme_classic() +
#   theme(text = element_text(size = 35))

# ggplot(out, aes(x = year, y = r2_ff6)) +
#   geom_line(lwd = 1) +
#   geom_hline(yintercept = 0, lty = 3) +
#   ggtitle("R2 FF6") +
#   theme_classic() +
#   theme(text = element_text(size = 35))



data <- copy(data_bk)

## Turnover dimension
# Calculate fund-level averages first, then merge back
fund_turnover <- data %>%
  group_by(wficn) %>%
  summarise(ave_turnover = mean(turnover, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    bin = ntile(ave_turnover, 5)
  ) %>%
  select(-ave_turnover) %>%
  setDT()

# Join back to original data
data <- data %>%
  left_join(fund_turnover,
    by = "wficn"
  )

# another way is to sort by period
data <- copy(data_bk)
data[, bin := ntile(turnover, 5), yyyymm]

fe <- fixef(sizefe)$wficn
data$phi <- rep(fe, as.numeric(table(data$wficn)))
mod1 <- lmer(alpha_CAPM ~ lag_logtna * as.factor(bin) + phi + si + (1 | wficn), data = data)
obj1 <- robust.se.remodel(data = data, model = mod1, cluster = c("wficn", "year"))
print(obj1)
