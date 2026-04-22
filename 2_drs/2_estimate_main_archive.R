# --- directly modified from Min's code, Analysis_DRS_forJ.R
library(this.path)
setwd(this.path::this.dir())
source("~/.runmefirst")
library(tidyverse)
library(lubridate)
library(lme4)
library(fixest)


########################## function definitions  ##########################################
source("../../../data_M/cluster2.R")

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


#######################################################
#### drs analysis
#######################################################

# load data
dat <- readRDS("../../../data/funds/min_sample.RDS")[order(fundid, yyyymm)] %>% setDT()
dat$year <- floor(dat$yyyymm / 100)
dat <- dat[order(fundid, yyyymm)]

## fill the missing ner (net expense ratio) with the latest ner information
dat <- dat |>
  group_by(fundid) |>
  fill(ner, .direction = "downup") |>
  ungroup() %>%
  setDT()

dat$ner <- DescTools::Winsorize(dat$ner, quantile(dat$ner, probs = c(0.01, 0.99)))

# mark if a fund is dead
max_ym <- dat[, max(yyyymm)]
dat <- dat |>
  group_by(fundid) |>
  mutate(dead = max(yyyymm) < max_ym, logage = log(fund_age), si = first(lag_logtna)) %>%
  setDT()
rm(max_ym)

### standardize q_{i0} by its category mean and variance
tmp <- dat |>
  group_by(fundid) |>
  summarise(q0 = first(lag_logtna), age0 = first(fund_age), B = first(morningstar_category))

category_mu_sd <- tmp |>
  group_by(B) |>
  summarise(q0mu.cat = mean(q0), q0sd.cat = sd(q0))

tmp <- tmp |>
  left_join(category_mu_sd, by = c("B")) |>
  mutate(zq0 = (q0 - q0mu.cat) / q0sd.cat) |>
  select(fundid, zq0)

dat <- dat |> left_join(tmp, by = c("fundid"))
dat_bk <- copy(dat)


## -------- size dynamics

size_models <- list()
size_models[[1]] <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner | fundid + yyyymm, cluster = c("fundid", "year"), data = dat)
size_models[[2]] <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + yyyymm, cluster = c("fundid", "year"), data = dat)
size_models[[3]] <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner | fundid, cluster = c("fundid", "year"), data = dat)
size_models[[4]] <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid, cluster = c("fundid", "year"), data = dat)



# save results
out <- data.table() # coefs
out_fe <- data.table() # fixed effects

for (i in 1:length(size_models)) {
  out <- rbind(out, data.table(
    spec_idx = i, var = names(size_models[[i]]$coefficients),
    coef = size_models[[i]]$coefficients, se = size_models[[i]]$se,
    obs = size_models[[i]]$nobs, ar2 = r2(size_models[[i]])["ar2"], war2 = r2(size_models[[i]])["war2"]
  ))

  fe_all <- fixef(size_models[[i]])
  phi_vector <- fe_all$fundid

  out_fe <- rbind(out_fe, data.table(
    spec_idx = i,
    fundid = names(phi_vector),
    phi = as.numeric(phi_vector)
  ))
}

to_dir <- "../tmp/drs/estimation/"
dir.create(to_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(out, paste0(to_dir, "size_dynamics.RDS"))
saveRDS(out_fe, paste0(to_dir, "size_dynamics_fe.RDS"))

## Table 4
################################################################################
### drs estimate by existing methods
p.drs <- function(data, ret.name = "benchmark_adj_gret") {
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

  # summarize and output
  out <- data.table(
    ret.name = ret.name,
    method = c("OLS", "FE", "RD", "RE"),
    coef = best,
    tstat = tv
  )

  return(out)
}


# estimate using various risk-adjusted dynamics (later parallelize)
ret_names <- c("benchmark_adj_gret", "capm_adj_gret", "ff3_adj_gret", "carhart_adj_gret", "ff5_adj_gret", "ff6_adj_gret")

tic()
plan(multisession, workers = parallel::detectCores() - 2)
out <- rbindlist(future_lapply(ret_names, function(x) p.drs(dat, x), future.seed = 123, future.packages = c("data.table")))
plan(sequential)
toc()

to_dir <- "../tmp/drs/estimation/"
dir.create(to_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(out, paste0(to_dir, "existing_methods.RDS"))


## ----- WB estimates
# femod <- size_models[[1]]
# i <- 1
# data <- dat
# ret.name <- 'benchmark_adj_gret'
p.wb <- function(size_models, i, dat = dat, ret.name = "benchmark_adj_gret") {
  femod <- size_models[[i]]
  int <- fixef(femod)$fundid
  dat$phi <- rep(int, as.numeric(table(dat$fundid)))
  setnames(dat, ret.name, "lhs")
  mod <- lmer(lhs ~ lag_logtna + phi + si + (1 | fundid), data = dat)
  r2_val <- suppressWarnings(performance::r2(mod))
  obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))

  out <- data.table(
    size_model = i, var = rownames(obj),
    ret.name = ret.name,
    coef = obj[, 1], se = obj[, 2], nobs = nobs(mod), marginalR2 = unlist(r2_val[2])
  )

  return(out)
}

vars <- c("benchmark_adj_gret", "capm_adj_gret", "ff3_adj_gret", "carhart_adj_gret", "ff5_adj_gret", "ff6_adj_gret")

out <- data.table()
for (this_var in vars) {
  tic(this_var)
  plan(multisession, workers = parallel::detectCores() - 2)
  out <- rbind(out, rbindlist(future_lapply(1:length(size_models), function(i) p.wb(size_models, i, dat, this_var), future.seed = 123, future.packages = c("data.table"))))
  plan(sequential)
  toc()
}

saveRDS(out, "../tmp/drs/estimation/wb.RDS")


# --- check variation of a_i and phi_i

dat <- copy(dat_bk)

# estimate phi and coefs
i <- 1
femod <- size_models[[i]]
int <- fixef(femod)$fundid
dat$phi <- rep(int, as.numeric(table(dat$fundid)))
mod <- lmer(benchmark_adj_gret ~ lag_logtna + phi + si + (1 | fundid), data = dat)
ss <- summary(mod)

dat[, ret_no_size := benchmark_adj_gret - ss$coefficients[2, 1] * lag_logtna]
tmp <- dat[, .(yyyymm, morningstar_category, fundid, ret_no_size, phi)]

# get rid of benchmark-based covariantion
tmp[, ret_no_size_demeaned := ret_no_size - mean(ret_no_size), yyyymm]
out <- tmp[, .(
  obs = .N, a = mean(ret_no_size_demeaned), se_a = sd(ret_no_size_demeaned) / sqrt(.N),
  phi = mean(phi), phi = last(phi)
), fundid]

# get summaries
p.summarize <- function(min_obs) {
  return(out[obs >= min_obs, .(
    min_obs = min_obs, num_funds = .N,
    sd_a_naive = sd(a),
    se_a = mean(se_a),
    sd_a_corrected = sqrt(sd(a)^2 - mean(se_a)^2),
    sd_phi = sd(phi), cor_a_phi = cor(a, phi)
  )])
}

options(width = 160)
rbindlist(lapply(c(0, 24, 60, 120, 180, 240, 300), p.summarize))

# what about the expected gross return at a point in time?
dat[, a := mean(ret_no_size), fundid]
dat[, expected_gross_return := a - ss$coefficients[2, 1] * lag_logtna]
dat[, sd(expected_gross_return), yyyymm][, mean(V1)]

dat[, cor(a, lag_logtna), yyyymm][, mean(V1)]
dat[, cor(a, lag_logtna, method = "spearman"), yyyymm][, mean(V1)]



# -- what if we use RD estimate?

dat <- copy(dat_bk)

i <- 1
femod <- size_models[[i]]
int <- fixef(femod)$fundid
dat$phi <- rep(int, as.numeric(table(dat$fundid)))

b <- -0.0026

dat[, ret_no_size := benchmark_adj_gret - b * lag_logtna]
tmp <- dat[, .(yyyymm, morningstar_category, fundid, ret_no_size, phi)]
tmp[, ret_no_size_demeaned := ret_no_size - mean(ret_no_size), yyyymm]
out <- tmp[, .(
  obs = .N, a = mean(ret_no_size_demeaned), se_a = sd(ret_no_size_demeaned) / sqrt(.N),
  phi = mean(phi), phi = last(phi)
), fundid]

rbindlist(lapply(c(0, 24, 60, 120, 180, 240, 300), p.summarize))

dat[, ret_no_size := benchmark_adj_gret + b * lag_logtna]

tmp <- dat[, .(yyyymm, fundid, ret_no_size)]
tmp[, ret_no_size_demeaned := ret_no_size - mean(ret_no_size), yyyymm]
out <- dat[, .(obs = .N, a = mean(ret_no_size_demeaned), sd_a = sd(ret_no_size_demeaned) / sqrt(.N), phi = mean(phi), phi = last(phi)), fundid]


# what about the expected gross return at a point in time?
dat[, a := mean(ret_no_size), fundid]
dat[, expected_gross_return := a - ss$coefficients[2, 1] * lag_logtna]
dat[, sd(expected_gross_return), yyyymm][, mean(V1)]

xx <- rnorm(1e5)
mean(xx > 1)
mean(xx < 1)




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
  # fmla1 = paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + ", perf),  "fundid+yyyymm", sep ="|")
  # fmla2 = paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd + ", perf),  "fundid+yyyymm", sep ="|")
  # fmla3 = paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + ", perf),  "fundid", sep ="|")
  # fmla4 = paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd + ", perf),  "fundid", sep ="|")
  if (sizemodel == 1) {
    fmla <- paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + ", perf), "fundid+yyyymm", sep = "|")
  } else if (sizemodel == 2) {
    fmla <- paste(paste0("logtna ~ lag_logtna + net_return + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd + ", perf), "fundid+yyyymm", sep = "|")
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
res <- wb.est.heter(data = dat1, sizefe = size_models[[1]])

## bin the data by 30%, 40%, 30%
dat1 <- singlesort(dat, nbins = NULL, ptiles = c(0.3, 0.7))
res <- wb.est.heter(data = dat1, sizefe = size_models[[1]], doublesort = F)

# doublesort
dat1 <- doublesort(dat, var1 = "activeness", var2 = "turnover", nbins = 3)
dat1 <- doublesort(dat, var1 = "turnover", var2 = "activeness", nbins = 3)
dat1 <- doublesort(dat, var1 = "ner", var2 = "activeness", nbins = 3)

res <- wb.est.heter(data = dat1, sizefe = size_models[[1]], doublesort = T)



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
  femod <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + yyyymm, cluster = c("fundid", "year"), data = dat)
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


sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + yyyymm, cluster = c("fundid", "year"), data = dat)
int <- fixef(sizefe1)$fundid
a <- c(summary(int), sd(int))

dd <- dat |>
  group_by(fundid) |>
  mutate(obs = n()) |>
  filter(obs >= 60) |>
  ungroup() |>
  dplyr::select(-obs)
sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + yyyymm, cluster = c("fundid", "year"), data = dd)
int <- fixef(sizefe1)$fundid
a6 <- c(summary(int), sd(int))

dd <- dat |>
  group_by(fundid) |>
  mutate(obs = n()) |>
  filter(obs >= 120) |>
  ungroup() |>
  dplyr::select(-obs)
sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + yyyymm, cluster = c("fundid", "year"), data = dd)
int <- fixef(sizefe1)$fundid
a12 <- c(summary(int), sd(int))

dd <- dat |>
  group_by(fundid) |>
  mutate(obs = n()) |>
  filter(obs >= 180) |>
  ungroup() |>
  dplyr::select(-obs)
sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + yyyymm, cluster = c("fundid", "year"), data = dd)
int <- fixef(sizefe1)$fundid
a18 <- c(summary(int), sd(int))

dd <- dat |>
  group_by(fundid) |>
  mutate(obs = n()) |>
  filter(obs >= 240) |>
  ungroup() |>
  dplyr::select(-obs)
sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + yyyymm, cluster = c("fundid", "year"), data = dd)
int <- fixef(sizefe1)$fundid
a24 <- c(summary(int), sd(int))

dd <- dat |>
  group_by(fundid) |>
  mutate(obs = n()) |>
  filter(obs >= 300) |>
  ungroup() |>
  dplyr::select(-obs)
sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + yyyymm, cluster = c("fundid", "year"), data = dd)
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
