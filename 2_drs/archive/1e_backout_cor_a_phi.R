# --- First get the point estimates of (a_i, phi_i) from the sample
########################################################################################
##  morningstar data analysis:  all analysis are based on gross port return
########################################################################################
library(this.path)
setwd(this.path::this.dir())
library(tidyverse)
library(lubridate)
library(lme4)
library(fixest)

# load data
load("../../../Data_M/MFdata.RData")

########################## function definitions  ##########################################
source("../../../Data_M/cluster2.R")

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


robust.se.remodel <- function(dat, model, cluster = c("fundid", "year"), weights = FALSE) {
  dd <- data.frame(y = model.frame(model)[, 1], model.matrix(model)[, -1])

  random_intercepts <- ranef(model)$fundid
  dd$raneff <- rep(random_intercepts[, 1], as.numeric(table(dat$fundid)))

  dd$y <- dd$y - dd$raneff
  xnames <- setdiff(names(dd)[-1], "raneff")
  fmla <- paste0("y ~ ", paste(xnames, collapse = " + "))
  # browser()
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

# robust.se.remodel(dat = dat, model=mod, cluster= c("fundid", "year"))
# robust.se.remodel(dat = dat, model=mod, cluster= c("morningstar_category", "year"))

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

category_mu_sd <- tmp |>
  group_by(B) |>
  summarise(q0mu.cat = mean(q0), q0sd.cat = sd(q0))

tmp <- tmp |>
  left_join(category_mu_sd, by = c("B")) |>
  mutate(zq0 = (q0 - q0mu.cat) / q0sd.cat) |>
  select(fundid, zq0)

dat <- dat |> left_join(tmp, by = c("fundid"))


### size process

sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner | fundid + caldt_end, cluster = c("fundid", "year"), data = dat)
sizefe3 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner | fundid, cluster = c("fundid", "year"), data = dat)

sizefe2 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + caldt_end, cluster = c("fundid", "year"), data = dat)
sizefe4 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid, cluster = c("fundid", "year"), data = dat)

## Table 2: size dynamics
summary(sizefe1)
summary(sizefe2)
summary(sizefe3)
summary(sizefe4)


## Table 4
################################################################################
### drs estimate by existing methods
drs <- function(data, ret.name = "benchmark_adj_gret") {
  dd <- data |> mutate(ret = !!sym(ret.name)) # convert ret.name to a symbol of data first (sym), then evaluate (!!) the symbol to a column:

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

o1 <- drs(data = dat, ret.name = "benchmark_adj_gret")
o2 <- drs(data = dat, ret.name = "capm_adj_gret")
o3 <- drs(data = dat, ret.name = "ff3_adj_gret")
o4 <- drs(data = dat, ret.name = "carhart_adj_gret")
o5 <- drs(data = dat, ret.name = "ff5_adj_gret")
o6 <- drs(data = dat, ret.name = "ff6_adj_gret")
o7 <- drs(data = dat, ret.name = "peer_adj_gret")

## drs by existing methods using different risk adjusted returns
t(rbind(o1$best, o2$best, o3$best, o4$best, o5$best, o6$best))
t(rbind(o1$tv, o2$tv, o3$tv, o4$tv, o5$tv, o6$tv))


## new method estimates
wb <- function(femod, data = dat, ni = 24) {
  int <- fixef(femod)$fundid
  dat$phi <- rep(int, as.numeric(table(dat$fundid)))
  mod <- lmer(benchmark_adj_gret ~ lag_logtna + phi + si + (1 | fundid), data = dat)
  obj1 <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
  # obj2 = robust.se.remodel(dat = dat, model=mod, cluster= c("fundid", "year"))
  return(obj1)
}

## CURRENT: trying to infer a and phi
# femod <- sizefe1
# data <- dat
back_out_a_and_phi <- function(femod, data = dat) {
  int <- fixef(femod)$fundid
  dat$phi <- rep(int, as.numeric(table(dat$fundid)))
  mod <- lmer(benchmark_adj_gret ~ lag_logtna + phi + si + (1 | fundid), data = dat)

  # back out a_i as average residual by fund
  dat$resid <- residuals(mod)
  data <- as.data.table(dat)

  out <- data[, .(obs = .N, a = mean(fixef(mod)[1] + fixef(mod)[3] * phi + resid), phi = last(phi)), by = fundid]
  out[, cor(a, phi)]
  out[, bin := ntile(obs, 20)]
  ggplot(out[, .(obs = mean(obs), cor = cor(a, phi)), by = bin], aes(x = obs, y = cor)) +
    geom_point(cex = 5) +
    theme(text = element_text(size = 35))

  obj1 <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
  # obj2 = robust.se.remodel(dat = dat, model=mod, cluster= c("fundid", "year"))
  return(obj1)
}


res1 <- wb(sizefe1)
res1[[1]]

res2 <- wb(sizefe2)
res3 <- wb(sizefe3)
res4 <- wb(sizefe4)

cbind(res1[, 1], res2[, 1], res3[, 1], res4[, 1]) # coef
cbind(res1[, 2], res2[, 2], res3[, 2], res4[, 2]) # tv

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

cbind(res10[, 1], res20[, 1], res30[, 1], res40[, 1]) # coef
cbind(res10[, 2], res20[, 2], res30[, 2], res40[, 2]) # tv


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

## function for subsample analysis
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


## function for full sample analysis (size X bins on various heterogeneity dimension)
wb.est.heter <- function(data, sizefe, doublesort = FALSE) {
  fe <- fixef(sizefe)$fundid
  # browser()
  # dat$phi[dat$fundid %in% names(fixef(fe)$fundid)] = rep(fe_intercepts, as.numeric(table(dat$fundid)))
  data$phi <- rep(fe, as.numeric(table(data$fundid)))

  if (doublesort) {
    mod <- lmer(benchmark_adj_gret ~ lag_logtna * double_bin + phi + si + (1 | fundid), data = data)
    obj <- robust.se.remodel(dat = data, model = mod, cluster = c("fundid", "year"))
    return(list(obj))
  } else {
    mod1 <- lmer(benchmark_adj_gret ~ lag_logtna * R2_bin + phi + si + (1 | fundid), data = data)
    mod2 <- lmer(benchmark_adj_gret ~ lag_logtna * ner_bin + phi + si + (1 | fundid), data = data)
    mod3 <- lmer(benchmark_adj_gret ~ lag_logtna * turnover_bin + phi + si + (1 | fundid), data = data)

    obj1 <- robust.se.remodel(dat = data, model = mod1, cluster = c("fundid", "year"))
    obj2 <- robust.se.remodel(dat = data, model = mod2, cluster = c("fundid", "year"))
    obj3 <- robust.se.remodel(dat = data, model = mod3, cluster = c("fundid", "year"))

    return(list(obj1, obj2, obj3))
  }
}


## prepare data: choose heterogeneity dimension and how many bins
# function to  bin values with equal or custom percentiles:
bin_by_percentile <- function(data, value_col, n_bins = NULL, percentiles = NULL,
                              bin_col_name = "bin") {
  # Check inputs
  if (is.null(n_bins) & is.null(percentiles)) {
    stop("Must provide either n_bins or percentiles")
  }

  if (!is.null(n_bins) & !is.null(percentiles)) {
    stop("Provide only one of n_bins or percentiles, not both")
  }

  # Create percentile breaks
  if (!is.null(n_bins)) {
    # Equal bins
    breaks <- seq(0, 1, length.out = n_bins + 1)
  } else {
    # Custom bins - percentiles should be cumulative (e.g., c(0.30, 0.70, 1.00))
    if (percentiles[1] != 0) {
      percentiles <- c(0, percentiles)
    }
    if (percentiles[length(percentiles)] != 1) {
      percentiles <- c(percentiles, 1)
    }
    breaks <- percentiles
  }

  # Calculate quantiles
  value_vec <- pull(data, {{ value_col }})
  quantile_breaks <- quantile(value_vec, probs = breaks, na.rm = TRUE)

  # Create bins
  # browser()
  data %>%
    mutate(!!sym(bin_col_name) := cut(value_vec,
      breaks = quantile_breaks,
      labels = 1:(length(breaks) - 1),
      include.lowest = TRUE
    ) %>% as.numeric())
}


# examples
# Equal bins (quintiles - 5 bins of 20% each)
# ave.turnover <- bin_by_percentile(ave.turnover, ave_turnover, n_bins = 5)
#
# # Custom bins (30%, 40%, 30%)
# ave.turnover <- bin_by_percentile(ave.turnover, ave_turnover,
#                                   percentiles = c(0.30, 0.70))
#
# # Custom bins with custom column name
# ave.turnover <- bin_by_percentile(ave.turnover, ave_turnover,
#                                   percentiles = c(0.30, 0.70),
#                                   bin_col_name = "turnover_tercile")

singlesort <- function(dat, nbins = NULL, ptiles = NULL) {
  ## Activeness dimension:  proxied by R2
  R2 <- dat |>
    group_by(fundid) |>
    summarise(activeness = sd(carhart_adj_gret) / sd(gross_return))
  # tmp <- tmp %>%
  #   mutate(R2_bin = ntile(activeness, 5))   ## 5 bins
  R2 <- bin_by_percentile(
    data = R2, value_col = "activeness", n_bins = nbins, percentiles = ptiles,
    bin_col_name = "R2_bin"
  )

  ## Turnover dimension
  ave.turnover <- dat |>
    group_by(fundid) |>
    summarise(ave_turnover = mean(turnover, na.rm = TRUE), .groups = "drop") |>
    filter(!is.na(ave_turnover))

  ave.turnover <- bin_by_percentile(
    data = ave.turnover, value_col = "ave_turnover", n_bins = nbins, percentiles = ptiles,
    bin_col_name = "turnover_bin"
  )


  ## Expense dimension
  ave.ner <- dat |>
    group_by(fundid) |>
    summarise(ave_ner = mean(ner, na.rm = TRUE), .groups = "drop")
  ave.ner <- bin_by_percentile(
    data = ave.ner, value_col = "ave_ner", n_bins = nbins, percentiles = ptiles,
    bin_col_name = "ner_bin"
  )

  dat |>
    left_join(R2 |> select(fundid, R2_bin), by = "fundid") |>
    left_join(ave.turnover |> select(fundid, turnover_bin), by = "fundid") |>
    left_join(ave.ner |> select(fundid, ner_bin), by = "fundid") |>
    mutate(
      turnover_bin = coalesce(turnover_bin, R2_bin), # Fill missing turnover with R2
      across(c(R2_bin, turnover_bin, ner_bin), as.factor)
    )

  return(dat)
}



## double sort, sort by var1 first, then with each var1_bin,  sort by var2
doublesort <- function(dat, var1, var2, nbins = 5) {
  # Create column names
  var1_bin_col <- paste0(var1, "_bin")
  var2_bin_col <- paste0(var2, "_bin")

  # First sort: bin by var1
  var1_summary <- dat |>
    group_by(fundid) |>
    summarise(value1 = case_when(
      var1 == "activeness" ~ sd(carhart_adj_gret) / sd(gross_return),
      var1 == "turnover" ~ mean(turnover, na.rm = TRUE),
      var1 == "ner" ~ mean(ner, na.rm = TRUE)
    ), .groups = "drop")

  # Fill missing turnover with sample mean
  if (var1 == "turnover") {
    mean_turnover <- mean(var1_summary$value1, na.rm = TRUE)
    var1_summary <- var1_summary |>
      mutate(value1 = if_else(is.na(value1), mean_turnover, value1))
  } else {
    var1_summary <- var1_summary |> filter(!is.na(value1))
  }

  var1_summary <- var1_summary |>
    bin_by_percentile(
      value_col = "value1", n_bins = nbins,
      bin_col_name = var1_bin_col
    )

  # Join first bin to data
  dat_with_bin1 <- dat |>
    left_join(var1_summary |> select(fundid, all_of(var1_bin_col)),
      by = "fundid"
    )

  # Second sort: within each bin of var1, bin by var2
  var2_summary <- dat_with_bin1 |>
    group_by(fundid, .data[[var1_bin_col]]) |>
    summarise(value2 = case_when(
      var2 == "activeness" ~ sd(carhart_adj_gret) / sd(gross_return),
      var2 == "turnover" ~ mean(turnover, na.rm = TRUE),
      var2 == "ner" ~ mean(ner, na.rm = TRUE)
    ), .groups = "drop")

  # Fill missing turnover with sample mean
  if (var2 == "turnover") {
    mean_turnover <- mean(var2_summary$value2, na.rm = TRUE)
    var2_summary <- var2_summary |>
      mutate(value2 = if_else(is.na(value2), mean_turnover, value2))
  } else {
    var2_summary <- var2_summary |> filter(!is.na(value2))
  }

  var2_summary <- var2_summary |>
    group_by(.data[[var1_bin_col]]) |>
    mutate(!!var2_bin_col := as.numeric(cut(
      value2,
      breaks = quantile(value2, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE),
      labels = 1:nbins,
      include.lowest = TRUE
    ))) |>
    ungroup()

  # Join both bins and create interaction bin
  dat <- dat |>
    left_join(var1_summary |> select(fundid, all_of(var1_bin_col)), by = "fundid") |>
    left_join(var2_summary |> select(fundid, all_of(var2_bin_col)), by = "fundid") |>
    mutate(
      double_bin = paste0(.data[[var1_bin_col]], "-", .data[[var2_bin_col]]),
      across(c(all_of(var1_bin_col), all_of(var2_bin_col), double_bin), as.factor)
    )

  # Set 1-1 as baseline - do this after ungrouping
  all_levels <- expand.grid(1:nbins, 1:nbins) |>
    arrange(Var1, Var2) |>
    mutate(level = paste0(Var1, "-", Var2)) |>
    pull(level)

  result <- dat |>
    mutate(double_bin = factor(double_bin, levels = all_levels))

  # levels = c("1-1", setdiff(sort(unique(double_bin)), "1-1"))))

  result
}


## single sort: bin the data by 5
dat1 <- singlesort(dat, nbins = 5, ptiles = NULL)
res <- wb.est.heter(data = dat1, sizefe = sizefe1)

## bin the data by 30%, 40%, 30%
dat1 <- singlesort(dat, nbins = NULL, ptiles = c(0.3, 0.7))
res <- wb.est.heter(data = dat1, sizefe = sizefe1, doublesort = F)



# doublesort
dat1 <- doublesort(dat, var1 = "activeness", var2 = "turnover", nbins = 3)
dat1 <- doublesort(dat, var1 = "turnover", var2 = "activeness", nbins = 3)
dat1 <- doublesort(dat, var1 = "ner", var2 = "activeness", nbins = 3)

res <- wb.est.heter(data = dat1, sizefe = sizefe1, doublesort = T)







## subsample

drs.highactive <- wb.est.heter.subsampe(dat[dat$R2_bin == 5, ], sizemodel = 1)
drs.lowactive <- wb.est.heter.subsampe(dat[dat$R2_bin == 1, ], sizemodel = 1)

drs.highactive
drs.lowactive

drs.highturnover <- wb.est.heter.subsampe(dat[dat$turnover_bin == 5, ], sizemodel = 1)
drs.lowturnover <- wb.est.heter.subsampe(dat[dat$turnover_bin == 1, ], sizemodel = 1)

drs.highturnover
drs.lowturnover

high <- drs(data = dat[dat$fundid %in% highto, ], ret.name = "benchmark_adj_gret")
low <- drs(dat[dat$fundid %in% lowto, ])

drs.highner <- wb.est.heter.subsampe(dat[dat$ner_bin == 5, ], sizemodel = 1)
drs.lowner <- wb.est.heter.subsampe(dat[dat$ner_bin == 1, ], sizemodel = 1)

drs.highner
drs.lowner

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


# ## Turnover dimension
# # Calculate fund-level averages first, then merge back
# fund_turnover <- dat %>%
#   group_by(fundid) %>%
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
#
# # Join back to original data
# dat <- dat %>%
#   left_join(fund_turnover %>% select(fundid, active),
#             by = "fundid")
#
#
# dat %>% count(liquidity, active)








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
