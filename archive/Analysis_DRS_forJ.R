########################################################################################
##  morningstar data analysis:  all analysis are based on gross port return
########################################################################################


# setwd("C:\\Users\\uqmzhu2\\Dropbox\\Jiacui Li\\diseconomy_of_scale\\Data_M")
library(this.path)
setwd(this.path::this.dir())
library(tidyverse)
library(lubridate)
library(lme4)
library(fixest)
library(dplyr)
library(data.table)

# load data
load("../../Data_M/MFdata.RData")

# # learn data
# data <- dat %>% as.data.table()
# tmp <- data[, .(obs = .N, tna = mean(tna, na.rm = TRUE)), caldt_end]
# tmp[, yyyy := year(caldt_end)]
# tmp <- tmp[order(caldt_end)][, .(obs = round(mean(obs)), tna = round(mean(tna))), yyyy]

########################## function definitions  ##########################################
source("../../Data_M/cluster2.R")

# function to estimate RD
estimation.rd <- function(dat) #  ZHU's method
{
  library(ivreg)
  library(sandwich)
  library(lmtest)
  rd <- ivreg(ret_RD ~ logx_RD + 0 | lag_logtna, data = dat)
  # Tests :   produce the same clustered standard errors as in Stata
  # b1 = coeftest.cluster(data=ddat, fm=rd, cluster1="SectorMonth", cluster2="FundId")      # double Clustered
  b1 <- coeftest.cluster(data = dat, fm = rd, cluster1 = "fundid", cluster2 = "year")

  # res = c(b1[, c(1,3,4)], rd$nobs)
  res <- c(b1[, c(1, 3, 4)])
  return(res)
}

# function to compute robust standard errors
robust.se.remodel <- function(dat, model, cluster = c("fundid", "year"), weights = FALSE) {
  dd <- data.frame(y = model.frame(model)[, 1], model.matrix(model)[, -1])

  random_intercepts <- ranef(model)$fundid
  dd$raneff <- rep(random_intercepts[, 1], as.numeric(table(dat$fundid)))

  dd$y <- dd$y - dd$raneff
  xnames <- setdiff(names(dd)[-1], "raneff")
  fmla <- paste0("y ~ ", paste(xnames, collapse = " + "))
  if (weights) {
    m1 <- lm(fmla, data = dd, weights = w)
  } else {
    m1 <- lm(fmla, data = dd)
  }

  dd$fundid <- dat$fundid
  dd$year <- dat$year
  # model.frame(model)  # Get the actual data used in fitting with the 1st column always the response variable
  # model.matrix(model)   # Get the model matrix (with interaction terms expanded)

  if (length(cluster) == 1) {
    return(coeftest.cluster.re(data = dd, fm = m1, cluster1 = cluster, weights = weights))
  }
  if (length(cluster) == 2) {
    return(coeftest.cluster.re(data = dd, fm = m1, cluster1 = cluster[1], cluster2 = cluster[2], weights = weights))
  }
}

# function to winsorize (J: I think this can just use DescTools::Winsorize)
winsorize <- function(x, q = 0.05, right = TRUE) {
  extrema <- quantile(x, c(q, 1 - q), na.rm = TRUE)
  if (right) {
    x[x > extrema[2]] <- extrema[2]
  } else {
    x[x < extrema[1]] <- extrema[1]
    x[x > extrema[2]] <- extrema[2]
  }
  x
}



#######################################################
#### drs analysis
#######################################################


## fill the missing ner (net expense ratio) with the latest ner information
dat <- dat |>
  group_by(fundid) |>
  fill(ner, .direction = "downup") |>
  ungroup()
dat$ner <- winsorize(dat$ner, q = 0.01, right = FALSE)

dat <- dat |>
  group_by(fundid) |>
  mutate(dead = max(caldt_end) < "2022-12-31", logage = log(fund_age), si = first(lag_logtna))

### standardize q_{i0} by its category mean and variance
tmp <- dat |>
  group_by(fundid) |>
  summarise(q0 = first(lag_logtna), age0 = first(fund_age), B = first(morningstar_category))
c(mean(tmp$q0), sd(tmp$q0), quantile(tmp$q0, c(1, 25, 50, 75, 99) / 100))
c(mean(tmp$age0), sd(tmp$age0), quantile(tmp$age0, c(1, 25, 50, 75, 99) / 100))
cor(tmp$q0, log(tmp$age0), method = "spearman")

# > summary(tmp$q0)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 2.708   3.046   3.794   4.258   5.140   9.884
# > summary(tmp$age0)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.08219  0.25205  0.91781  3.74537  3.25205 56.50685

# variation of iniitial size by category
category_mu_sd <- tmp |>
  group_by(B) |>
  summarise(q0mu.cat = mean(q0), q0sd.cat = sd(q0))

# turn initial size into z-score
tmp <- tmp |>
  left_join(category_mu_sd, by = c("B")) |>
  mutate(zq0 = (q0 - q0mu.cat) / q0sd.cat) |>
  select(fundid, zq0)

dat <- dat |> left_join(tmp, by = c("fundid"))
#
#
# # dat = dat |> mutate(across(gross_return, ~ map_dfc(1:12, ~ lag(.gross_return, .y)) |> set_names(paste0("gret_lag", 1:12)))) |> ungroup()  # error in code
# dat = dat |>  mutate(
#   x = benchmark_adj_gret,  #net_return,   # benchmark_ret,
#   x_lag1 = lag(x, 1),
#   x_lag2 = lag(x, 2),
#   x_lag3 = lag(x, 3),
#   x_lag4 = lag(x, 4),
#   x_lag5 = lag(x, 5),
#   x_lag6 = lag(x, 6),
#   x_lag7 = lag(x, 7),
#   x_lag8 = lag(x, 8),
#   x_lag9 = lag(x, 9),
#   x_lag10 = lag(x, 10),
#   x_lag11 = lag(x, 11),
#   x_lag12 = lag(x, 12)
# )  |> ungroup()


### UNDERSTAND FLOW PROCESS
# dat = dat |> group_by(fundid) |> mutate(flow1 = logtna - lag_logtna - net_return,
#                                         flow2 = tna_cpi/exp(lag_logtna)-(1+net_return))
# dat$flow = dat$flow1
# dat$flow = winsorize(dat$flow2, q=0.01, right=FALSE)
#
# p = 1/100
# quantile(dat$flow, c(p, 1-p), na.rm=T)
#
# flowfe = feols(flow ~ lag_logtna + net_return + expanding_ret + expanding_ret_sd + expanding_bret + expanding_bret_sd + logage + ner
#                  | fundid + caldt_end, cluster = c("fundid", "year"), data=dat)
#
# flowlm = lm(flow ~ lag_logtna + net_return + expanding_ret + expanding_ret_sd + expanding_bret + expanding_bret_sd + logage + ner, data =dat)
# sizelm1 = lm(logtna ~ lag_logtna + net_return + expanding_ret + expanding_ret_sd + expanding_bret + expanding_bret_sd + logage + ner, data=dat)
# sizelm2 = lm(logtna ~ lag_logtna + net_return + expanding_ret + expanding_bret, data=dat)
#
#
# int <- fixef(flowfe)$fundid
# int <- fixef(flowfe)$caldt_end
# c(summary(int), sd(int))


### estimate size process (different fixed effect choices)

sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner | fundid + caldt_end, cluster = c("fundid", "year"), data = dat)
sizefe3 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner | fundid, cluster = c("fundid", "year"), data = dat)

sizefe2 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + caldt_end, cluster = c("fundid", "year"), data = dat)
sizefe4 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid, cluster = c("fundid", "year"), data = dat)

## Table 2: size dynamics
summary(sizefe1)
summary(sizefe2)
summary(sizefe3)
summary(sizefe4)

# ### this size process is too complicated. many lagged alpha measures are not very useful.
# sizefe3 = feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd
#                + x_lag1+ x_lag2 + x_lag3+ x_lag4+ x_lag5+ x_lag6+ x_lag7+ x_lag8+ x_lag9+ x_lag10+ x_lag10+ x_lag12
#                | fundid+caldt_end, cluster = c("fundid", "year"), data=dat)
# sizefe31= feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd
#                 + x_lag1+ x_lag2 + x_lag3+ x_lag4+ x_lag5+ x_lag6+ x_lag7+ x_lag8+ x_lag9+ x_lag10+ x_lag10+ x_lag12
#                 | fundid, cluster = c("fundid", "year"), data=dat)



## Table 4
################################################################################
### drs estimate by existing methods
drs <- function(data, ret.name = "benchmark_adj_gret") { # browser()
  dd <- data |> mutate(ret = !!sym(ret.name)) # convert ret.name to a symbol of data first (sym), then evaluate (!!) the symbol to a column:

  # RD
  dd_rd <- dd |>
    group_by(fundid) |>
    mutate(
      ret_RD = ret - rev(cumsum(rev(ret))) / c(n():1),
      logx_RD = lag_logtna - rev(cumsum(rev(lag_logtna))) / c(n():1)
    ) |>
    filter(!(ret_RD == 0 & logx_RD == 0)) |>
    ungroup()

  ## model fitting
  ols <- lm(ret ~ lag_logtna, data = dd) # ols
  obj1 <- coeftest.cluster(data = dd, fm = ols, cluster1 = c("fundid"), cluster2 = c("year"))

  fe <- feols(ret ~ lag_logtna | fundid, cluster = c("fundid", "year"), data = dd) ## fe

  rdout <- estimation.rd(dat = dd_rd) # rd

  re <- lmer(ret ~ lag_logtna + (1 | fundid), data = dd) # re
  obj2 <- robust.se.remodel(dat = dd, model = re, cluster = c("fundid", "year"))

  best <- c(coef(ols)[2], coef(fe), rdout[1], fixef(re)[2])
  tv <- c(obj1[2, 3], fe$coeftable[, 3], rdout[2], obj2[2, 3])

  return(list(best = best, tv = tv))
  # return(list(ols=ols, fe=fe, rd = rdout, re=re))
}

# different ways to adjust returns
o1 <- drs(data = dat, ret.name = "benchmark_adj_gret")
o2 <- drs(data = dat, ret.name = "capm_adj_gret")
o3 <- drs(data = dat, ret.name = "ff3_adj_gret")
o4 <- drs(data = dat, ret.name = "carhart_adj_gret")
o5 <- drs(data = dat, ret.name = "ff5_adj_gret")
o6 <- drs(data = dat, ret.name = "ff6_adj_gret")
o7 <- drs(data = dat, ret.name = "peer_adj_gret")

## drs by existing methods using different risk adjusted returns (so there is some variation)
# rows: ols, fe, rd, re
round(t(rbind(o1$best, o2$best, o3$best, o4$best, o5$best, o6$best)) * 1e4, 2)
round(t(rbind(o1$tv, o2$tv, o3$tv, o4$tv, o5$tv, o6$tv)), 2)

## new method estimates
wb <- function(femod, data = dat, ni = 24) {
  int <- fixef(femod)$fundid
  dat$phi <- rep(int, as.numeric(table(dat$fundid)))
  mod <- lmer(benchmark_adj_gret ~ lag_logtna + phi + si + (1 | fundid), data = dat)
  # return( fixef(mod)[2])
  # VarCorr(mod)
  # browser()
  obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
  return(obj)
}

# different size processes
res1 <- wb(sizefe1)
res2 <- wb(sizefe2)
res3 <- wb(sizefe3)
res4 <- wb(sizefe4)

# round(cbind(res1[[2]][2, 1], res2[[2]][2, 1], res3[[2]][2, 1], res4[[2]][2, 1]) * 1e4, 2)
# round(cbind(res1[[2]][2, 2], res2[[2]][2, 2], res3[[2]][2, 1], res4[[2]][2, 1]) * 1e4, 2)

round(cbind(res1[, 1], res2[, 1], res3[, 1], res4[, 1]) * 1e4, 2) # coef
round(cbind(res1[, 3], res2[, 3], res3[, 3], res4[, 3]), 2) # tv - not as high power as we have hoped

## wb on standardized initial size (standerdised by style)
wb.zq0 <- function(femod) {
  int <- fixef(femod)$fundid
  dat$phi <- rep(int, as.numeric(table(dat$fundid)))
  mod <- lmer(benchmark_adj_gret ~ lag_logtna + phi + zq0 + (1 | fundid), data = dat)
  # return( fixef(mod)[2])
  obj1 <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
  return(obj1[, -c(2, 4)])
}
res10 <- wb.zq0(sizefe1)
res20 <- wb.zq0(sizefe2)
res30 <- wb.zq0(sizefe3)
res40 <- wb.zq0(sizefe4)

round(cbind(res10[, 1], res20[, 1], res30[, 1], res40[, 1]) * 1e4, 2) # coef
round(cbind(res10[, 2], res20[, 2], res30[, 2], res40[, 2]), 2) # tv


wb_drop_q0 <- function(femod) {
  int <- fixef(femod)$fundid
  dat$phi <- rep(int, as.numeric(table(dat$fundid)))
  mod <- lmer(benchmark_adj_gret ~ lag_logtna + phi + (1 | fundid), data = dat)
  # return( fixef(mod)[2])
  obj1 <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
  return(obj1[, -c(2, 4)])
}

wb_drop_q0(sizefe1)



######################################################################################################################
####### heterogenous drs
######################################################################################################################

## subsample analysis
wb.est.heter.subsampe <- function(data, ni = 12, ret.name = "benchmark_adj_gret", sizemodel = 1) {
  dat <- data |>
    group_by(fundid) |>
    mutate(obs = n()) |>
    filter(obs >= ni) |>
    ungroup() |>
    dplyr::select(-obs)
  dat <- dat |> mutate(ret = !!sym(ret.name))

  if (ret.name == "benchmark_adj_gret") {
    perf <- "expanding_ret"
  } else if (ret.name == "capm_adj_gret") {
    perf <- "expanding_capm"
  } else if (ret.name == "ff3_adj_gret") {
    perf <- "expanding_ff3"
  } else if (ret.name == "carhart_adj_gret") {
    perf <- "expanding_car"
  } else if (ret.name == "ff5_adj_gret") {
    perf <- "expanding_ff5"
  } else {
    perf <- "expanding_ff6"
  }

  # four size process corresponding to the ones reported in the paper
  # fmla1 = paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + ", perf),  "fundid+caldt_end", sep ="|")
  # fmla2 = paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd + ", perf),  "fundid+caldt_end", sep ="|")
  # fmla3 = paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + ", perf),  "fundid", sep ="|")
  # fmla4 = paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd + ", perf),  "fundid", sep ="|")

  if (sizemodel == 1) {
    fmla <- paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + ", perf), "fundid+caldt_end", sep = "|")
  } else if (sizemodel == 2) {
    fmla <- paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd + ", perf), "fundid+caldt_end", sep = "|")
  } else if (sizemodel == 3) {
    fmla <- paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + ", perf), "fundid", sep = "|")
  } else {
    fmla <- paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd + ", perf), "fundid", sep = "|")
  }

  fe <- feols(as.formula(fmla), cluster = c("fundid", "year"), data = dat)
  fe_intercepts1 <- fixef(fe)$fundid

  # dat$phi[dat$fundid %in% names(fixef(fe)$fundid)] = rep(fe_intercepts, as.numeric(table(dat$fundid)))
  dat$phi <- rep(fe_intercepts1, as.numeric(table(dat$fundid)))
  mod1 <- lmer(ret ~ lag_logtna + phi + si + (1 | fundid), data = dat)
  # mod1 = lmer(ret ~ I(exp(lag_logtna)) + phi + si + (1 | fundid), data=dat)
  obj1 <- robust.se.remodel(dat = dat, model = mod1, cluster = c("fundid", "year"))

  # return(list(obj1,fe))
  return(obj1[, -c(2, 4)])
}

## Activeness dimension:  proxied by R2
tmp <- dat |>
  group_by(fundid) |>
  summarise(activeness = sd(carhart_adj_gret) / sd(gross_return))
cutoff <- median(tmp$activeness)
highactive <- tmp$fundid[tmp$activeness >= cutoff]
lowactive <- tmp$fundid[tmp$activeness < cutoff]

drs.highactive <- wb.est.heter.subsampe(dat[dat$fundid %in% highactive, ], sizemodel = 1)
drs.lowactive <- wb.est.heter.subsampe(dat[dat$fundid %in% lowactive, ], sizemodel = 1)

drs.highactive <- wb.est.heter.subsampe(dat[dat$fundid %in% highactive, ], sizemodel = 4)
drs.lowactive <- wb.est.heter.subsampe(dat[dat$fundid %in% lowactive, ], sizemodel = 4)

drs.highactive
drs.lowactive

high <- drs(data = dat[dat$fundid %in% highactive, ], ret.name = "benchmark_adj_gret")
low <- drs(dat[dat$fundid %in% lowactive, ])

## Turnover dimension
dat <- dat |>
  group_by(fundid) |>
  mutate(ave_turnover = mean(turnover, na.rm = TRUE))
tmp <- dat[!is.na(dat$ave_turnover), ]
ave.turnover <- tmp |>
  group_by(fundid) |>
  summarise(ave_turnover = first(ave_turnover))
cutoff <- median(ave.turnover$ave_turnover)
highto <- tmp$fundid[tmp$ave_turnover >= cutoff]
lowto <- tmp$fundid[tmp$ave_turnover < cutoff]

drs.highturnover <- wb.est.heter.subsampe(dat[dat$fundid %in% highto, ], sizemodel = 1)
drs.lowturnover <- wb.est.heter.subsampe(dat[dat$fundid %in% lowto, ], sizemodel = 1)

drs.highturnover
drs.lowturnover

high <- drs(data = dat[dat$fundid %in% highto, ], ret.name = "benchmark_adj_gret")
low <- drs(dat[dat$fundid %in% lowto, ])


## Liquidity dimension
# new estimator
drs.large <- wb.est.heter.subsampe(dat[dat$morningstar_category %in% c("US Fund Large Growth", "US Fund Large Blend", "US Fund Large Value"), ])
drs.mid <- wb.est.heter.subsampe(dat[dat$morningstar_category %in% c("US Fund Mid-Cap Growth", "US Fund Mid-Cap Blend", "US Fund Mid-Cap Value"), ])
drs.small <- wb.est.heter.subsampe(dat[dat$morningstar_category %in% c("US Fund Small Growth", "US Fund Small Blend", "US Fund Small Value"), ])

# existing estimate
large <- drs(data = dat[dat$morningstar_category %in% c("US Fund Large Growth", "US Fund Large Blend", "US Fund Large Value"), ], ret.name = "benchmark_adj_gret")
mid <- drs(dat[dat$morningstar_category %in% c("US Fund Mid-Cap Growth", "US Fund Mid-Cap Blend", "US Fund Mid-Cap Value"), ])
small <- drs(dat[dat$morningstar_category %in% c("US Fund Small Growth", "US Fund Small Blend", "US Fund Small Value"), ])



######################### fullsample, two-way interaction of size and portfolio characteristics

dat <- dat %>%
  mutate(
    liquidity = case_when(
      morningstar_category %in% c("US Fund Large Growth", "US Fund Large Blend", "US Fund Large Value") ~ "large",
      morningstar_category %in% c("US Fund Mid-Cap Growth", "US Fund Mid-Cap Blend", "US Fund Mid-Cap Value") ~ "mid",
      morningstar_category %in% c("US Fund Small Growth", "US Fund Small Blend", "US Fund Small Value") ~ "small",
      TRUE ~ NA_character_ # for any other categories not specified
    ),
    liquidity = factor(liquidity, levels = c("large", "mid", "small"))
  )

dat %>% count(morningstar_category, liquidity)


## Turnover dimension
# Calculate fund-level averages first, then merge back
fund_turnover <- dat %>%
  group_by(fundid) %>%
  summarise(ave_turnover = mean(turnover, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    median_turnover = median(ave_turnover, na.rm = TRUE),
    active = if_else(
      ave_turnover >= median_turnover,
      "high",
      "low"
    ),
    active = factor(active, levels = c("low", "high"))
  )

# Join back to original data
dat <- dat %>%
  left_join(fund_turnover %>% select(fundid, active),
    by = "fundid"
  )


dat %>% count(liquidity, active)


wb.est.heter <- function(data, sizefe) {
  fe <- fixef(sizefe)$fundid

  # dat$phi[dat$fundid %in% names(fixef(fe)$fundid)] = rep(fe_intercepts, as.numeric(table(dat$fundid)))
  data$phi <- rep(fe, as.numeric(table(data$fundid)))
  mod1 <- lmer(benchmark_adj_gret ~ lag_logtna * liquidity + phi + si + (1 | fundid), data = data)

  data2 <- data %>% filter(!is.na(active)) # there is missing active info
  mod2 <- lmer(benchmark_adj_gret ~ lag_logtna * active + phi + si + (1 | fundid), data = data2)
  # mod1 = lmer(benchmark_adj_gret ~ lag_logtna*liquidity*active - liquiditymid:activehigh - liquiditysmall:activehigh+ phi + si + (1 | fundid), data=data)
  # mod1 = lmer(benchmark_adj_gret ~ liquidity + active + lag_logtna + lag_logtna:liquidity+ lag_logtna:active + lag_logtna:liquidity:active+ phi + si + (1 | fundid), data=data)
  # browser()
  obj1 <- robust.se.remodel(dat = data, model = mod1, cluster = c("fundid", "year"))
  obj2 <- robust.se.remodel(dat = data2, model = mod2, cluster = c("fundid", "year"))

  # return(list(obj1,obj2))
  return(list(mod1, mod2))
}

res <- wb.est.heter(data = dat, sizefe = sizefe1)

# FE
feols(benchmark_adj_gret ~ lag_logtna * liquidity | fundid, dat, cluster = c("fundid", "year"))
feols(benchmark_adj_gret ~ lag_logtna * active | fundid, dat, cluster = c("fundid", "year"))

# RE
m1 <- lmer(benchmark_adj_gret ~ lag_logtna * liquidity + (1 | fundid), dat)
robust.se.remodel(dat = dat, model = m1, cluster = c("fundid", "year"))

dat2 <- dat %>% filter(!is.na(active)) # there is missing active info
m2 <- lmer(benchmark_adj_gret ~ lag_logtna * active + (1 | fundid), dat2)
robust.se.remodel(dat = dat2, model = m2, cluster = c("fundid", "year"))


## subsample analysis
drs(data = dat[!is.na(dat$active) & dat$active == "high", ])











### estimate drs with longer time series
wb.est(data = dat, ret.name = "benchmark_adj_gret", ni = 120)

wb.est <- function(dat) {
  femod <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + caldt_end, cluster = c("fundid", "year"), data = dat)
  int <- fixef(femod)$fundid
  dat$phi <- rep(int, as.numeric(table(dat$fundid)))
  mod <- lmer(benchmark_adj_gret ~ lag_logtna + phi + si + (1 | fundid), data = dat)
  # mod = lmer(benchmark_adj_gret ~ lag_logtna + phi + (1 | fundid), data=dat)
  # return( fixef(mod)[2])
  obj1 <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
  return(obj1[, -c(2, 4)])
}

dd <- dat |>
  filter(dead == 0) |>
  group_by(fundid) |>
  mutate(obs = n()) |>
  filter(obs >= 24) |>
  ungroup() |>
  dplyr::select(-obs) ## survivor funds
dd <- dat |>
  group_by(fundid) |>
  mutate(obs = n()) |>
  filter(obs >= 60) |>
  ungroup() |>
  dplyr::select(-obs)
length(unique(dd$fundid))
dim(dd)[1]
# dd = dat |> group_by(fundid) |> mutate(obs = n()) |> filter(obs >= 120, fund_age > 6) |> ungroup()  |> dplyr::select(-obs)
# drs(dat,ret.name="benchmark_adj_gret")
drs(dd, ret.name = "benchmark_adj_gret")
wb.est(dd)


sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + caldt_end, cluster = c("fundid", "year"), data = dat)
int <- fixef(sizefe1)$fundid
a <- c(summary(int), sd(int))

dd <- dat |>
  group_by(fundid) |>
  mutate(obs = n()) |>
  filter(obs >= 60) |>
  ungroup() |>
  dplyr::select(-obs)
sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + caldt_end, cluster = c("fundid", "year"), data = dd)
int <- fixef(sizefe1)$fundid
a6 <- c(summary(int), sd(int))

dd <- dat |>
  group_by(fundid) |>
  mutate(obs = n()) |>
  filter(obs >= 120) |>
  ungroup() |>
  dplyr::select(-obs)
sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + caldt_end, cluster = c("fundid", "year"), data = dd)
int <- fixef(sizefe1)$fundid
a12 <- c(summary(int), sd(int))

dd <- dat |>
  group_by(fundid) |>
  mutate(obs = n()) |>
  filter(obs >= 180) |>
  ungroup() |>
  dplyr::select(-obs)
sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + caldt_end, cluster = c("fundid", "year"), data = dd)
int <- fixef(sizefe1)$fundid
a18 <- c(summary(int), sd(int))

dd <- dat |>
  group_by(fundid) |>
  mutate(obs = n()) |>
  filter(obs >= 240) |>
  ungroup() |>
  dplyr::select(-obs)
sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + caldt_end, cluster = c("fundid", "year"), data = dd)
int <- fixef(sizefe1)$fundid
a24 <- c(summary(int), sd(int))

dd <- dat |>
  group_by(fundid) |>
  mutate(obs = n()) |>
  filter(obs >= 300) |>
  ungroup() |>
  dplyr::select(-obs)
sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + caldt_end, cluster = c("fundid", "year"), data = dd)
int <- fixef(sizefe1)$fundid
a30 <- c(summary(int), sd(int))

rbind(a, a6, a12, a18, a24, a30)


## check correlation between size fixed effects and a_i
mod1 <- wb(sizefe2)

a2 <- dat |>
  group_by(fundid) |>
  mutate(ait = benchmark_adj_gret - mod1[2, 1] * lag_logtna) |>
  summarise(a = mean(ait), s = first(si))
int <- fixef(sizefe2)$fundid
cor(int, a2$a)
cor(int, a2$s)
cor(a2$s, a2$a)




fe1 <- feols(logtna ~ lag_logtna + gross_return + expanding_ret + expanding_bret | fundid, cluster = c("fundid", "year"), data = dat)
int <- fixef(fe2)$fundid
dat$phi <- rep(int, as.numeric(table(dat$fundid)))
wb4 <- lmer(benchmark_adj_gret ~ lag_logtna + phi + (1 | fundid), data = dat)
summary(wb4)

cor(dat$phi, dat$age, method = "spearman")


re <- lmer(benchmark_adj_gret ~ lag_logtna + (1 | fundid), data = dat)
fe <- feols(benchmark_adj_gret ~ lag_logtna | fundid, cluster = c("fundid", "year"), data = dat)
fe_intercepts <- fixef(fe)$fundid
plot(density(fe_intercepts)) # skew-normal distribution
shapiro.test(fe_intercepts)

dat$phi <- rep(fe_intercepts, as.numeric(table(dat$fundid)))
wb4 <- lmer(benchmark_adj_gret ~ lag_logtna + phi + (1 | fundid), data = dat) # failed.  cannot fit a distribution

re <- lmer(benchmark_adj_gret ~ lag_logtna + (1 | fundid), data = dat)
re1 <- lmer(benchmark_adj_gret ~ lag_logtna + si + (1 | fundid), data = dat)
re2 <- lmer(benchmark_adj_gret ~ lag_logtna + expanding_bret + expanding_bret_sd + ner + (1 | fundid), data = dat)



## fund by fund ols
b <- dat |>
  group_by(fundid) |>
  summarise(
    b = coef(lm(capm_adj_gret ~ lag_logtna, data = pick(everything())))[2],
    a = coef(lm(capm_adj_gret ~ lag_logtna, data = pick(everything())))[1],
    sigma = summary(lm(capm_adj_gret ~ lag_logtna, data = pick(everything())))$sigma,
    n = n()
  )

weighted.mean(b$b, w = 1 / b$sigma^2)
weighted.mean(b$b, w = 1 / b$n)
