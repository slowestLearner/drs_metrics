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

########################## function definitions  ##########################################

# load some base functions (for clustering)
source("../../../data_M/cluster2.R")

# function to estimate RD
estimation.rd <- function(data) #  ZHU's method
{
  # data <- dd_rd
  library(ivreg)
  library(sandwich)
  library(lmtest)
  rd <- ivreg(ret_RD ~ logx_RD + 0 | lag_logtna, data = data)
  # Tests :   produce the same clustered standard errors as in Stata
  b1 <- coeftest.cluster(data = data, fm = rd, cluster1 = "wficn", cluster2 = "yyyymm")

  # res = c(b1[, c(1,3,4)], rd$nobs)
  res <- c(b1[, c(1, 3, 4)])
  return(res)
}

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
#### drs analysis
#######################################################

# data: merge with fund size, tc
# data <- readRDS("../tmp/raw_data/fund_net_ret_risk_adjusted.RDS")
data <- readRDS("../tmp/raw_data/fund_net_ret_risk_adjusted_full_sample.RDS")
tmp <- readRDS("../../../data/funds/active_equity_fund_monthly_data.RDS")[, .(yyyymm, wficn, exp_ratio, lag_logtna = log(mtna_1), turnover = turn_ratio)]
tmp[, idx := frank(yyyymm, ties.method = "dense")]
tmp <- tmp[tmp[, .(idx = idx - 1, wficn, logtna = lag_logtna)], on = .(idx, wficn), nomatch = NULL][, idx := NULL]
data <- data[tmp, on = .(yyyymm, wficn), nomatch = NULL]
data <- data[exp_ratio < .05] # otherwise a bit weird
rm(tmp)

# sometimes cluster by year instead of year-month (seems to make little difference)
data[, yyyy := floor(yyyymm / 100)]

# infer fund age from sample start
tmp <- unique(data[, .(yyyymm)])
tmp[, date := as.Date(as.character(yyyymm * 100 + 28), "%Y%m%d")]
data <- data[tmp, on = .(yyyymm)]
data[, min_date := min(date), by = wficn]
data[, fund_age := as.numeric(date - min_date) / 365][, min_date := NULL]
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

## Table 4
################################################################################
### drs estimate by existing methods
drs <- function(data, ret.name = "benchmark_adj_gret", time_cluster = "yyyymm") { # browser()
  dd <- data |> mutate(ret = !!sym(ret.name)) # convert ret.name to a symbol of data first (sym), then evaluate (!!) the symbol to a column:

  # data for RD
  dd_rd <- dd |>
    group_by(wficn) |>
    mutate(
      ret_RD = ret - rev(cumsum(rev(ret))) / c(n():1),
      logx_RD = lag_logtna - rev(cumsum(rev(lag_logtna))) / c(n():1)
    ) |>
    filter(!(ret_RD == 0 & logx_RD == 0)) |>
    ungroup()

  ## model fitting
  ols <- lm(ret ~ lag_logtna, data = dd) # ols
  obj1 <- coeftest.cluster(data = dd, fm = ols, cluster1 = c("wficn"), cluster2 = time_cluster)

  fe <- feols(ret ~ lag_logtna | wficn, cluster = c("wficn", "yyyymm"), data = dd) ## fe

  rdout <- estimation.rd(data = dd_rd) # rd
  # browser()
  dd <- dd |>
    group_by(wficn) |>
    mutate(obs = n()) |>
    filter(obs >= 24) |>
    ungroup() |>
    dplyr::select(-obs)
  re <- lmer(ret ~ lag_logtna + (1 | wficn), data = dd) # re
  obj2 <- robust.se.remodel(data = dd, model = re, cluster = c("wficn", time_cluster))

  best <- c(coef(ols)[2], coef(fe), rdout[1], fixef(re)[2])
  tv <- c(obj1[2, 3], fe$coeftable[, 3], rdout[2], obj2[2, 3])

  return(list(best = best, tv = tv))
}


# estimate DRS, with different ways to adjust returns
data[, year := floor(yyyymm / 100)]

tic()
o1 <- drs(data = data, ret.name = "alpha_CAPM")
o2 <- drs(data = data, ret.name = "alpha_FF3")
o3 <- drs(data = data, ret.name = "alpha_FFC")
o4 <- drs(data = data, ret.name = "alpha_FF5")
o5 <- drs(data = data, ret.name = "alpha_FF6")
gc()
toc()

## drs by existing methods using different risk adjusted returns (so there is some variation)
# rows: ols, fe, rd, re
round(t(rbind(o1$best, o2$best, o3$best, o4$best, o5$best)) * 1e2, 2)
round(t(rbind(o1$tv, o2$tv, o3$tv, o4$tv, o5$tv)), 2)

## within-between estimation - TO MIN: I keep getting "isSingular" warning in "lmer"
wb <- function(femod, data = data, ret.name = "alpha_CAPM", time_cluster = "yyyymm") {
  int <- fixef(femod)$wficn
  data$phi <- rep(int, as.numeric(table(data$wficn)))
  mod <- lmer(as.formula(paste0(ret.name, " ~ lag_logtna + phi + si + (1 | wficn)")), data = data)
  obj <- robust.se.remodel(data = data, model = mod, cluster = c("wficn", time_cluster))
  return(obj)
}

# different size processes - FOR MIN - I get isSingular warning
wb(sizefe1, data, "alpha_CAPM")
wb(sizefe1, data, "alpha_FF3")




######################################################################################################################
####### heterogenous drs
######################################################################################################################

## subsample analysis
wb.est.heter.subsampe <- function(data, ni = 12, ret.name = "alpha_CAPM", sizemodel = 1) {
  data <- data |>
    group_by(wficn) |>
    mutate(obs = n()) |>
    filter(obs >= ni) |>
    ungroup() |>
    dplyr::select(-obs)
  data <- data |> mutate(ret = !!sym(ret.name))

  # if (ret.name == "benchmark_adj_gret") {
  #   perf <- "expanding_ret"
  # } else if (ret.name == "capm_adj_gret") {
  #   perf <- "expanding_capm"
  # } else if (ret.name == "ff3_adj_gret") {
  #   perf <- "expanding_ff3"
  # } else if (ret.name == "carhart_adj_gret") {
  #   perf <- "expanding_car"
  # } else if (ret.name == "ff5_adj_gret") {
  #   perf <- "expanding_ff5"
  # } else {
  #   perf <- "expanding_ff6"
  # }

  # four size process corresponding to the ones reported in the paper
  if (sizemodel == 1) {
    fmla <- paste(paste0("logtna ~ lag_logtna + ret + logage + exp_ratio"), "wficn+year", sep = "|")
  } else if (sizemodel == 2) {
    fmla <- paste(paste0("logtna ~ lag_logtna + ret + logage + exp_ratio + expanding_ret_sd + expanding_bret + expanding_bret_sd"), "wficn+year", sep = "|")
  } else if (sizemodel == 3) {
    fmla <- paste(paste0("logtna ~ lag_logtna + ret + logage + exp_ratio"), "wficn", sep = "|")
  } else {
    fmla <- paste(paste0("logtna ~ lag_logtna + ret + logage + exp_ratio + expanding_ret_sd + expanding_bret + expanding_bret_sd"), "wficn", sep = "|")
  }

  fe <- feols(as.formula(fmla), cluster = c("wficn", "year"), data = data)
  fe_intercepts1 <- fixef(fe)$wficn


  # data$phi[data$wficn %in% names(fixef(fe)$wficn)] = rep(fe_intercepts, as.numeric(table(data$wficn)))
  data$phi <- rep(fe_intercepts1, as.numeric(table(data$wficn)))
  mod1 <- lmer(ret ~ lag_logtna + phi + si + (1 | wficn), data = data)
  obj1 <- robust.se.remodel(data = data, model = mod1, cluster = c("wficn", "year"))

  # return(list(obj1,fe))
  return(obj1[, -c(2, 4)])
}

## Activeness dimension:  proxied by R2
data <- copy(data_bk)
tmp <- data |>
  group_by(wficn) |>
  summarise(activeness = sd(alpha_FFC) / sd(ret)) %>%
  na.omit()
cutoff <- median(tmp$activeness)
highactive <- tmp$wficn[tmp$activeness >= cutoff]
lowactive <- tmp$wficn[tmp$activeness < cutoff]

# this is hardcoded. currently, wb.est.heter.subsampe clusters by year instead of year-month
data[, year := floor(yyyymm / 100)]

# estimate separately
this_ret_var <- "alpha_CAPM"

tic()
drs.highactive <- wb.est.heter.subsampe(data[data$wficn %in% highactive, ], ret.name = this_ret_var, sizemodel = 1)
drs.lowactive <- wb.est.heter.subsampe(data[data$wficn %in% lowactive, ], ret.name = this_ret_var, sizemodel = 1)
toc()

# print results
print(drs.highactive)
print(drs.lowactive)

# unpaired t-test
coef_diff <- drs.highactive[2, 1] - drs.lowactive[2, 1]
se.highactive <- drs.highactive[2, 1] / drs.highactive[2, 2]
se.lowactive <- drs.lowactive[2, 1] / drs.lowactive[2, 2]
coef_diff / sqrt(se.highactive^2 + se.lowactive^2)


## Turnover dimension
data <- data |>
  group_by(wficn) |>
  mutate(ave_turnover = mean(turnover, na.rm = TRUE))
tmp <- data[!is.na(data$ave_turnover), ]
ave.turnover <- tmp |>
  group_by(wficn) |>
  summarise(ave_turnover = first(ave_turnover))
cutoff <- median(ave.turnover$ave_turnover)
highto <- tmp$wficn[tmp$ave_turnover >= cutoff]
lowto <- tmp$wficn[tmp$ave_turnover < cutoff]

this_ret_var <- "alpha_CAPM"
drs.highturnover <- wb.est.heter.subsampe(data[data$wficn %in% highto, ], ret.name = this_ret_var, sizemodel = 1)
drs.lowturnover <- wb.est.heter.subsampe(data[data$wficn %in% lowto, ], ret.name = this_ret_var, sizemodel = 1)

drs.highturnover
drs.lowturnover

# unpaired t-test
coef_diff <- drs.highturnover[2, 1] - drs.lowturnover[2, 1]
se.highturnover <- drs.highturnover[2, 1] / drs.highturnover[2, 2]
se.lowturnover <- drs.lowturnover[2, 1] / drs.lowturnover[2, 2]
coef_diff / sqrt(se.highturnover^2 + se.lowturnover^2)

# ## Liquidity dimension
# # new estimator
# drs.large <- wb.est.heter.subsampe(data[data$morningstar_category %in% c("US Fund Large Growth", "US Fund Large Blend", "US Fund Large Value"), ])
# drs.mid <- wb.est.heter.subsampe(data[data$morningstar_category %in% c("US Fund Mid-Cap Growth", "US Fund Mid-Cap Blend", "US Fund Mid-Cap Value"), ])
# drs.small <- wb.est.heter.subsampe(data[data$morningstar_category %in% c("US Fund Small Growth", "US Fund Small Blend", "US Fund Small Value"), ])

# # existing estimate
# large <- drs(data = data[data$morningstar_category %in% c("US Fund Large Growth", "US Fund Large Blend", "US Fund Large Value"), ], ret.name = "benchmark_adj_gret")
# mid <- drs(data[data$morningstar_category %in% c("US Fund Mid-Cap Growth", "US Fund Mid-Cap Blend", "US Fund Mid-Cap Value"), ])
# small <- drs(data[data$morningstar_category %in% c("US Fund Small Growth", "US Fund Small Blend", "US Fund Small Value"), ])



# ######################### fullsample, two-way interaction of size and portfolio characteristics

# data <- data %>%
#   mutate(
#     liquidity = case_when(
#       morningstar_category %in% c("US Fund Large Growth", "US Fund Large Blend", "US Fund Large Value") ~ "large",
#       morningstar_category %in% c("US Fund Mid-Cap Growth", "US Fund Mid-Cap Blend", "US Fund Mid-Cap Value") ~ "mid",
#       morningstar_category %in% c("US Fund Small Growth", "US Fund Small Blend", "US Fund Small Value") ~ "small",
#       TRUE ~ NA_character_ # for any other categories not specified
#     ),
#     liquidity = factor(liquidity, levels = c("large", "mid", "small"))
#   )

# data %>% count(morningstar_category, liquidity)


# ## Turnover dimension
# # Calculate fund-level averages first, then merge back
# fund_turnover <- data %>%
#   group_by(wficn) %>%
#   summarise(ave_turnover = mean(turnover, na.rm = TRUE)) %>%
#   ungroup() %>%
#   mutate(
#     median_turnover = median(ave_turnover, na.rm = TRUE),
#     active = if_else(
#       ave_turnover >= median_turnover,
#       "high",
#       "low"
#     ),
#     active = factor(active, levels = c("low", "high"))
#   )

# # Join back to original data
# data <- data %>%
#   left_join(fund_turnover %>% select(wficn, active),
#     by = "wficn"
#   )


# data %>% count(liquidity, active)


# wb.est.heter <- function(data, sizefe) {
#   fe <- fixef(sizefe)$wficn

#   # data$phi[data$wficn %in% names(fixef(fe)$wficn)] = rep(fe_intercepts, as.numeric(table(data$wficn)))
#   data$phi <- rep(fe, as.numeric(table(data$wficn)))
#   mod1 <- lmer(benchmark_adj_gret ~ lag_logtna * liquidity + phi + si + (1 | wficn), data = data)

#   data2 <- data %>% filter(!is.na(active)) # there is missing active info
#   mod2 <- lmer(benchmark_adj_gret ~ lag_logtna * active + phi + si + (1 | wficn), data = data2)
#   # mod1 = lmer(benchmark_adj_gret ~ lag_logtna*liquidity*active - liquiditymid:activehigh - liquiditysmall:activehigh+ phi + si + (1 | wficn), data=data)
#   # mod1 = lmer(benchmark_adj_gret ~ liquidity + active + lag_logtna + lag_logtna:liquidity+ lag_logtna:active + lag_logtna:liquidity:active+ phi + si + (1 | wficn), data=data)
#   # browser()
#   obj1 <- robust.se.remodel(data = data, model = mod1, cluster = c("wficn", "year"))
#   obj2 <- robust.se.remodel(data = data2, model = mod2, cluster = c("wficn", "year"))

#   # return(list(obj1,obj2))
#   return(list(mod1, mod2))
# }

# res <- wb.est.heter(data = data, sizefe = sizefe1)

# # FE
# feols(benchmark_adj_gret ~ lag_logtna * liquidity | wficn, data, cluster = c("wficn", "year"))
# feols(benchmark_adj_gret ~ lag_logtna * active | wficn, data, cluster = c("wficn", "year"))

# # RE
# m1 <- lmer(benchmark_adj_gret ~ lag_logtna * liquidity + (1 | wficn), data)
# robust.se.remodel(data = data, model = m1, cluster = c("wficn", "year"))

# data2 <- data %>% filter(!is.na(active)) # there is missing active info
# m2 <- lmer(benchmark_adj_gret ~ lag_logtna * active + (1 | wficn), data2)
# robust.se.remodel(data = data2, model = m2, cluster = c("wficn", "year"))


# ## subsample analysis
# drs(data = data[!is.na(data$active) & data$active == "high", ])











# ### estimate drs with longer time series
# wb.est(data = data, ret.name = "benchmark_adj_gret", ni = 120)

# wb.est <- function(data) {
#   femod <- feols(logtna ~ lag_logtna + ret + logage + exp_ratio + expanding_ret_sd + expanding_bret + expanding_bret_sd | wficn + yyyymm, cluster = c("wficn", "year"), data = data)
#   int <- fixef(femod)$wficn
#   data$phi <- rep(int, as.numeric(table(data$wficn)))
#   mod <- lmer(benchmark_adj_gret ~ lag_logtna + phi + si + (1 | wficn), data = data)
#   # mod = lmer(benchmark_adj_gret ~ lag_logtna + phi + (1 | wficn), data=data)
#   # return( fixef(mod)[2])
#   obj1 <- robust.se.remodel(data = data, model = mod, cluster = c("wficn", "year"))
#   return(obj1[, -c(2, 4)])
# }

# dd <- data |>
#   filter(dead == 0) |>
#   group_by(wficn) |>
#   mutate(obs = n()) |>
#   filter(obs >= 24) |>
#   ungroup() |>
#   dplyr::select(-obs) ## survivor funds
# dd <- data |>
#   group_by(wficn) |>
#   mutate(obs = n()) |>
#   filter(obs >= 60) |>
#   ungroup() |>
#   dplyr::select(-obs)
# length(unique(dd$wficn))
# dim(dd)[1]
# # dd = data |> group_by(wficn) |> mutate(obs = n()) |> filter(obs >= 120, fund_age > 6) |> ungroup()  |> dplyr::select(-obs)
# # drs(data,ret.name="benchmark_adj_gret")
# drs(dd, ret.name = "benchmark_adj_gret")
# wb.est(dd)


# sizefe1 <- feols(logtna ~ lag_logtna + ret + logage + exp_ratio + expanding_ret_sd + expanding_bret + expanding_bret_sd | wficn + yyyymm, cluster = c("wficn", "year"), data = data)
# int <- fixef(sizefe1)$wficn
# a <- c(summary(int), sd(int))

# dd <- data |>
#   group_by(wficn) |>
#   mutate(obs = n()) |>
#   filter(obs >= 60) |>
#   ungroup() |>
#   dplyr::select(-obs)
# sizefe1 <- feols(logtna ~ lag_logtna + ret + logage + exp_ratio + expanding_ret_sd + expanding_bret + expanding_bret_sd | wficn + yyyymm, cluster = c("wficn", "year"), data = dd)
# int <- fixef(sizefe1)$wficn
# a6 <- c(summary(int), sd(int))

# dd <- data |>
#   group_by(wficn) |>
#   mutate(obs = n()) |>
#   filter(obs >= 120) |>
#   ungroup() |>
#   dplyr::select(-obs)
# sizefe1 <- feols(logtna ~ lag_logtna + ret + logage + exp_ratio + expanding_ret_sd + expanding_bret + expanding_bret_sd | wficn + yyyymm, cluster = c("wficn", "year"), data = dd)
# int <- fixef(sizefe1)$wficn
# a12 <- c(summary(int), sd(int))

# dd <- data |>
#   group_by(wficn) |>
#   mutate(obs = n()) |>
#   filter(obs >= 180) |>
#   ungroup() |>
#   dplyr::select(-obs)
# sizefe1 <- feols(logtna ~ lag_logtna + ret + logage + exp_ratio + expanding_ret_sd + expanding_bret + expanding_bret_sd | wficn + yyyymm, cluster = c("wficn", "year"), data = dd)
# int <- fixef(sizefe1)$wficn
# a18 <- c(summary(int), sd(int))

# dd <- data |>
#   group_by(wficn) |>
#   mutate(obs = n()) |>
#   filter(obs >= 240) |>
#   ungroup() |>
#   dplyr::select(-obs)
# sizefe1 <- feols(logtna ~ lag_logtna + ret + logage + exp_ratio + expanding_ret_sd + expanding_bret + expanding_bret_sd | wficn + yyyymm, cluster = c("wficn", "year"), data = dd)
# int <- fixef(sizefe1)$wficn
# a24 <- c(summary(int), sd(int))

# dd <- data |>
#   group_by(wficn) |>
#   mutate(obs = n()) |>
#   filter(obs >= 300) |>
#   ungroup() |>
#   dplyr::select(-obs)
# sizefe1 <- feols(logtna ~ lag_logtna + ret + logage + exp_ratio + expanding_ret_sd + expanding_bret + expanding_bret_sd | wficn + yyyymm, cluster = c("wficn", "year"), data = dd)
# int <- fixef(sizefe1)$wficn
# a30 <- c(summary(int), sd(int))

# rbind(a, a6, a12, a18, a24, a30)


# ## check correlation between size fixed effects and a_i
# mod1 <- wb(sizefe2)

# a2 <- data |>
#   group_by(wficn) |>
#   mutate(ait = benchmark_adj_gret - mod1[2, 1] * lag_logtna) |>
#   summarise(a = mean(ait), s = first(si))
# int <- fixef(sizefe2)$wficn
# cor(int, a2$a)
# cor(int, a2$s)
# cor(a2$s, a2$a)




# fe1 <- feols(logtna ~ lag_logtna + gross_return + expanding_bret | wficn, cluster = c("wficn", "year"), data = data)
# int <- fixef(fe2)$wficn
# data$phi <- rep(int, as.numeric(table(data$wficn)))
# wb4 <- lmer(benchmark_adj_gret ~ lag_logtna + phi + (1 | wficn), data = data)
# summary(wb4)

# cor(data$phi, data$age, method = "spearman")


# re <- lmer(benchmark_adj_gret ~ lag_logtna + (1 | wficn), data = data)
# fe <- feols(benchmark_adj_gret ~ lag_logtna | wficn, cluster = c("wficn", "year"), data = data)
# fe_intercepts <- fixef(fe)$wficn
# plot(density(fe_intercepts)) # skew-normal distribution
# shapiro.test(fe_intercepts)

# data$phi <- rep(fe_intercepts, as.numeric(table(data$wficn)))
# wb4 <- lmer(benchmark_adj_gret ~ lag_logtna + phi + (1 | wficn), data = data) # failed.  cannot fit a distribution

# re <- lmer(benchmark_adj_gret ~ lag_logtna + (1 | wficn), data = data)
# re1 <- lmer(benchmark_adj_gret ~ lag_logtna + si + (1 | wficn), data = data)
# re2 <- lmer(benchmark_adj_gret ~ lag_logtna + expanding_bret + expanding_bret_sd + exp_ratio + (1 | wficn), data = data)



# ## fund by fund ols
# b <- data |>
#   group_by(wficn) |>
#   summarise(
#     b = coef(lm(capm_adj_gret ~ lag_logtna, data = pick(everything())))[2],
#     a = coef(lm(capm_adj_gret ~ lag_logtna, data = pick(everything())))[1],
#     sigma = summary(lm(capm_adj_gret ~ lag_logtna, data = pick(everything())))$sigma,
#     n = n()
#   )

# weighted.mean(b$b, w = 1 / b$sigma^2)
# weighted.mean(b$b, w = 1 / b$n)
