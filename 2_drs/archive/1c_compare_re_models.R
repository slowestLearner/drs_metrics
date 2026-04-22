# --- updated Min's code to use nlme instead. Sanity check that results are the same
library(this.path)
setwd(this.path::this.dir())
library(tidyverse)
library(lubridate)
library(lme4)
library(fixest)
library(dplyr)
library(data.table)
library(nlme)

# load data
load("../../../Data_M/MFdata.RData")

# source("../../Data_M/cluster2.R")

# function to compute robust standard errors: updated
robust.se.remodel_j <- function(dat, model, cluster = c("fundid", "year"), weights = FALSE) {
  # --- BRIDGE FOR NLME VS LMER ---
  if (inherits(model, "lme")) {
    # 1. Reconstruct y (Response)
    y_val <- model$residuals[, "fixed"] + model$fitted[, "fixed"]

    # 2. Extract X (Design Matrix)
    # We get the fixed formula (e.g., y ~ x1 + x2)
    fixed_formula <- formula(model)

    # We use model.frame to handle NAs correctly based on the model's own data
    m_frame <- model.frame(fixed_formula, data = model$data)
    X_mat <- model.matrix(fixed_formula, data = m_frame)

    # Remove the intercept column (usually the 1st) because your
    # robust.se.remodel_j logic builds its own or handles it in lm()
    X_mat <- X_mat[, -1, drop = FALSE]

    random_intercepts <- nlme::ranef(model)
  } else if (inherits(model, "lmerMod")) {
    # Original lmer logic
    y_val <- model.frame(model)[, 1]
    X_mat <- model.matrix(model)[, -1]
    random_intercepts <- lme4::ranef(model)$fundid
  }

  # Create the data frame for the second-stage OLS
  dd <- data.frame(y = y_val, X_mat)

  # Handle Random Intercepts
  # table(dat$fundid) works if dat is sorted and matches the model frame
  dd$raneff <- rep(random_intercepts[, 1], as.numeric(table(dat$fundid)))

  # De-mean the response variable by the random effects
  dd$y <- dd$y - dd$raneff

  xnames <- setdiff(names(dd), c("y", "raneff"))
  fmla <- as.formula(paste0("y ~ ", paste(xnames, collapse = " + ")))

  # Fit the OLS model
  if (weights) {
    # Assumes 'w' exists in 'dat' or 'dd'
    m1 <- lm(fmla, data = dd, weights = dat$w)
  } else {
    m1 <- lm(fmla, data = dd)
  }

  dd$fundid <- dat$fundid
  dd$year <- dat$year

  # Standard error clustering
  if (length(cluster) == 1) {
    return(coeftest.cluster.re(data = dd, fm = m1, cluster1 = cluster[1], weights = weights))
  }
  if (length(cluster) == 2) {
    return(coeftest.cluster.re(data = dd, fm = m1, cluster1 = cluster[1], cluster2 = cluster[2], weights = weights))
  }
}

# multiply
dat$benchmark_adj_gret <- dat$benchmark_adj_gret * 100

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

### estimate size process (different fixed effect choices)
sizefe1 <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner | fundid + caldt_end, cluster = c("fundid", "year"), data = dat)



# new code
wb_multi_method <- function(femod, dat = dat, method = "nlme") {
  int <- fixef(femod)$fundid
  dat$phi <- rep(int, as.numeric(table(dat$fundid)))

  # Ensure fundid is a factor for grouping
  data <- as.data.table(dat)
  data[, fundid := as.factor(fundid)]

  # Define the formula string
  form_str <- "benchmark_adj_gret ~ lag_logtna + phi + si"

  if (method == "lmer") {
    # Original lme4 approach
    mod <- lmer(as.formula(paste0(form_str, " + (1 | fundid)")), data = data)

    # Using your existing robust SE helper
    rr <- robust.se.remodel_j(dat = data, model = mod, cluster = c("fundid", "year"))
  } else if (method == "nlme") {
    # nlme approach - slower but allows for AR(1) or heteroskedasticity
    mod <- lme(as.formula(form_str), data = data, random = ~ 1 | fundid)
    # rr <- robust.se.remodel_j(dat = data, model = mod, cluster = c("fundid", "year"))
    # rr <- summary(mod)
    rr <- robust.se.remodel_j(dat = data, model = mod, cluster = c("fundid", "year"))
  }
  return(data.table(method, var = rownames(rr), coef = rr[, 1], se = rr[, 2], tstat = rr[, 3]))
}

out <- rbindlist(lapply(c("lmer", "nlme"), function(method) wb_multi_method(sizefe1, dat, method = method)))
out[var == "lag_logtna"]
