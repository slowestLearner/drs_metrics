# --- changed nothing, except to parallelize it and output in csv. Goal is to roughly replicate earlier results
library(this.path)
setwd(this.path::this.dir())
source("../../../Data_M/cluster2.R")
library(tictoc)
library(tidyverse)
library(fixest)
library(data.table)
# setFixest_notes(FALSE)
library(lme4)
library(MASS)
library(tictoc)
library(future.apply)

# gendata is for the scenario of both growth (a_i, phi_i) and inception(a_i, q_{i0}) matching.
gendata <- function(gm = 0.5, im = 0, b = 0, N = 1000, tt = 200, phi.mean, phi.sd, cor_a_stde = 0.15) { ### b: drs
  ######### Size equation: size_{it} = phi_i + rho1*size_{it-1} + rho2*(B_t + alpha_{it}) + rho3*alpha_{it-1} + rho4*B_{t-1} + rnorm(N, 0, stdv);

  n_years <- ceiling(tt / 12)

  # coefs in the size model
  rho1 <- 0.99 # 0.986
  rho2 <- 1.06 # 1.061
  rho3 <- 2.30 # 2.177
  rho4 <- 0 # 0.738
  stdq <- 0.09 # 0.2

  ######### Return equation:  alpha_{it} = a_i + b*size_{it-1} + epsilon_{it}
  ## 1) distribution  of a (parameters for monthly fund skill)
  amu <- 0.004 # 0.004
  asd <- 0.0036 # 0.0015   ## control the cross-section skill heterogeneity
  a <- rnorm(N, amu, asd) # monthly a value
  # a = runif(N, -0.001, 0.008)   # for robustness check:  when RE condition is not met

  ## 2) residual distribution: epsilon_{it}
  # stde = 0.02
  stde <- rgamma(N, shape = 4.1, rate = 225.6)

  ## 3)  large alpha is paired with large idiosyncratic volatility
  # cor(a, stde)
  # specify a correlation matrix
  # 0.15: target correlation for cor(a, stde)
  corr_mat <- matrix(cor_a_stde, ncol = 2, nrow = 2)

  diag(corr_mat) <- 1
  mvdat <- mvrnorm(n = N, mu = c(0, 0), Sigma = corr_mat, empirical = TRUE) # from library(MASS)
  rx <- rank(mvdat[, 1], ties.method = "first")
  ry <- rank(mvdat[, 2], ties.method = "first")
  a_sorted <- sort(a)
  stde_sorted <- sort(stde)
  # use the ranks to specify the order of the sorted data
  ae <- data.frame(a = a_sorted[rx], stde = stde_sorted[ry])
  # a = a_sorted[rx]
  # stde = stde_sorted[ry]
  # cor(ae$a, ae$stde) # the obtained correlation does not perfectly match the specified one, but the difference is relatively small.

  ######### Benchmark return distribution
  Bmu <- 0.01
  Bsd <- 0.05
  nstyle <- 5 # the number of styles / benchmarks
  Bret1 <- matrix(rnorm(nstyle * tt, Bmu, Bsd), nrow = nstyle, ncol = tt)
  cc <- N / nstyle
  Bret <- do.call(rbind, replicate(cc, Bret1, simplify = FALSE)) # matrix(rep(Bret1, each = cc), ncol = ncol(Bret1), byrow = FALSE)

  expanding_ave <- t(apply(Bret1[, -tt], 1, FUN = function(y) map_dbl(seq_along(y), ~ mean(y[1:.x])))) # dimension is nstyle X (tt-1)
  expanding_bret_ave <- do.call(rbind, replicate(cc, expanding_ave, simplify = FALSE)) # dimension is N X (tt-1)

  ## fund size grow rates phi_i
  phi <- rnorm(N, phi.mean, sd = phi.sd) # fund-specific growth effects for size
  # gm: target growth matching correlation
  corr_mat <- matrix(gm, ncol = 2, nrow = 2)
  diag(corr_mat) <- 1
  mvdat <- mvrnorm(n = N, mu = c(0, 0), Sigma = corr_mat, empirical = TRUE) # from library(MASS)
  rx <- rank(mvdat[, 1], ties.method = "first")
  ry <- rank(mvdat[, 2], ties.method = "first")
  # use the ranks to specify the order of the sorted phi while preserve the correlation of a and stde
  ae <- ae[order(ae$a), ]
  ae <- ae[rx, ]
  phi_sorted <- sort(phi)
  ae$phi <- phi_sorted[ry]

  ## Initial fund size
  # q0 = ifelse(samestartsize, log(rep(200,N)),  rnorm(N, log(200), 1) )  # Initial fund size in logarithm
  q0 <- runif(N, log(20), log(100)) # Initial fund size in logarithm

  # im: target correlation between q0 and a_i
  if (im != 0) {
    corr_mat <- matrix(im, ncol = 2, nrow = 2)
    diag(corr_mat) <- 1
    mvdat <- mvrnorm(n = N, mu = c(0, 0), Sigma = corr_mat, empirical = TRUE) # from library(MASS)
    rx <- rank(mvdat[, 1], ties.method = "first")
    ry <- rank(mvdat[, 2], ties.method = "first")

    ae <- ae[order(ae$a), ]
    ae <- ae[rx, ]
    q0_sorted <- sort(q0)
    q0 <- q0_sorted[ry]
  }

  a <- ae$a
  stde <- ae$stde
  phi <- ae$phi

  # print(c(cor(a, q0), cor(a, phi), cor(a, stde)))

  r <- q <- ps <- matrix(NA, nrow = N, ncol = tt) # ps:  perceived alpha/skill
  for (i in 1:tt) {
    if (i == 1) {
      r[, i] <- a + b * q0 + rnorm(N, 0, stde)
      q[, i] <- phi + rho1 * q0 + rho2 * (Bret[, i] + r[, i]) + rnorm(N, 0, stdq)
      ps[, i] <- r[, i]
    } else {
      r[, i] <- a + b * q[, i - 1] + rnorm(N, 0, stde)
      ps[, i] <- rowMeans(r[, c(1:i)])
      q[, i] <- phi + rho1 * q[, i - 1] + rho2 * (Bret[, i] + r[, i]) + rho3 * ps[, i - 1] + rho4 * expanding_bret_ave[, i - 1] + rnorm(N, 0, stdq)
      q[, i] <- pmax(q[, i], 1)
    }
  }

  # expanding mean for r: represent skill at each given t based on data upto t-1.
  perceived_r <- cbind(NA, ps[, -tt])

  ## 10 year cumulative survival rates taken from Linnainmaa (2013)
  sr.c <- c(0.995, 0.976, 0.943, 0.908, 0.868, 0.83, 0.789, 0.754, 0.717, 0.68)

  # if num of years change, then will have to extrapolate, etc


  ## yearly survival rate
  sr <- c(sr.c[1], sr.c[-1] / sr.c[-10])
  if (n_years > 10) {
    sr <- c(sr, rep(sr[length(sr)], n_years - 10))
  }
  if (n_years < 10) {
    sr <- sr[1:n_years]
  }




  ## alpha-based fund exit
  # tmp = fundsurvial(Nfunds=N, maxT = tt, alpha=r, sigma = stde, size=q)
  # idx1 = tmp$idx   # stopping time: a vector of length N
  # survrate1 = tmp$survrate
  # idx1[is.na(idx1)] = tt

  ## exiting rule based on alpha
  ### Bayesian updating of fund alpha
  m <- matrix(rep(NA, N * (tt + 1)), nrow = N)
  v <- matrix(rep(NA, N * (tt + 1)), nrow = N)
  sigma2 <- stde^2 ## monthly fund alpha variance

  m[, 1] <- 0 # alpha          ## prior:  mean of alpha dist
  v[, 1] <- (0.0125 / 12)^2 ## prior:  monthly var of alpha dist based on prior of annual alpha sd being 0.0125
  for (i in 2:(tt + 1)) {
    w <- v[, i - 1] / (v[, i - 1] + sigma2)
    m[, i] <- (1 - w) * m[, i - 1] + w * r[, i - 1]
    v[, i] <- sigma2 * w
  }
  m <- m[, -1] # posterior alpha
  v <- v[, -1] # posterior alpha variance

  idx1 <- rep(tt, N) # total survival months for each fund
  idx2 <- rep(tt, N)
  idx3 <- rep(tt, N)
  livefunds.alpha <- 1:N # alpha-based fund exit
  livefunds.size <- 1:N # size-based fund exit
  livefunds.random <- 1:N # random fund exit
  for (yr in 1:n_years) {
    ii <- ((yr - 1) * 12 + 1):((yr) * 12)

    ## posterior alpha (m)-based
    val <- rowMeans(m[livefunds.alpha, ii])
    quantile_values <- quantile(val, probs = 1 - sr[yr])
    dead <- livefunds.alpha[val <= quantile_values]
    idx1[dead] <- apply(m[dead, ii], 1, which.min) + (yr - 1) * 12
    livefunds.alpha <- setdiff(livefunds.alpha, dead)

    ## relative size(change)-based
    Q <- cbind(q0, q)
    delta_q <- Q[, -1] - Q[, -ncol(Q)] # monthly percentage change in fund TNA
    val <- rowMeans(delta_q[livefunds.size, ii])
    quantile_values <- quantile(val, probs = 1 - sr[yr])
    deadfunds <- livefunds.size[val <= quantile_values]
    idx2[deadfunds] <- apply(q[deadfunds, ii], 1, which.min) + (yr - 1) * 12
    livefunds.size <- setdiff(livefunds.size, deadfunds)

    ## random exit
    deadfunds.random <- sample(livefunds.random, length(deadfunds))
    idx3[deadfunds.random] <- sample(1:12, length(deadfunds), replace = TRUE) + (yr - 1) * 12
    livefunds.random <- setdiff(livefunds.random, deadfunds.random)
  }



  qq <- cbind(q0, q[, -tt]) # lagged log fund size
  # dd = data.frame(fundid= rep(1:N, each=tt), ret=as.vector(t(r)), logx=as.vector(t(qq)), survival =as.vector(t(survival)) )
  dd <- data.frame(
    fundid = rep(1:N, each = tt), ret = as.vector(t(r)), Bret = as.vector(t(Bret)), ps = as.vector(t(perceived_r)),
    bret_ave_exp = as.vector(t(cbind(NA, expanding_bret_ave))),
    a = rep(a, each = tt), phi = rep(phi, each = tt),
    lagsize = as.vector(t(qq)), size = as.vector(t(q)),
    idx1 = rep(idx1, each = tt), idx2 = rep(idx2, each = tt), idx3 = rep(idx3, each = tt)
  )
  # dd = data.frame(fundid= rep(1:N, each=tt), ret=as.vector(t(r)), logx=as.vector(t(qq)), idx = rep(idx, each=tt) )
  dd <- dd |>
    group_by(fundid) |>
    mutate(
      survival1 = ifelse(row_number() <= first(idx1), 1, 0),
      survival2 = ifelse(row_number() <= first(idx2), 1, 0),
      survival3 = ifelse(row_number() <= first(idx3), 1, 0)
    ) |>
    ungroup() |>
    dplyr::select(-idx1, -idx2, -idx3)

  dd$obs <- rep(c(1:tt), N)

  return(list(data = dd))
}


# TODO: replace with a standardized function
estimation.rd <- function(dat, weights = FALSE) #  ZHU's method
{
  library(ivreg)
  library(sandwich)
  library(lmtest)
  if (weights) {
    rd <- ivreg(ret_RD ~ logx_RD + 0 | lagsize, data = dat, weights = w)
  } else {
    rd <- ivreg(ret_RD ~ logx_RD + 0 | lagsize, data = dat)
  }
  # Tests :   produce the same clustered standard errors as in Stata
  # b1 = coeftest.cluster(data=ddat, fm=rd, cluster1="SectorMonth", cluster2="FundId")      # double Clustered
  b1 <- coeftest.cluster(data = dat, fm = rd, cluster1 = "fundid")

  # res = c(b1[, c(1,3,4)], rd$nobs)
  res <- c(b1[, c(1, 3, 4)])
  return(res)
}

# TODO: replace with a standardized function
robust.se.remodel <- function(dat, model, cluster = c("fundid")) {
  random_intercepts <- ranef(model)$fundid
  dat$raneff <- rep(random_intercepts[, 1], as.numeric(table(dat$fundid)))

  dat$y <- dat$ret - dat$raneff
  xnames <- names(fixef(model))[-1]
  fmla <- paste0("y ~ ", paste(xnames, collapse = " + "))
  mod1 <- lm(fmla, data = dat)
  if (length(cluster) == 1) {
    return(coeftest.cluster.re(data = dat, fm = mod1, cluster1 = cluster))
  }
  if (length(cluster) == 2) {
    return(coeftest.cluster.re(data = dat, fm = mod1, cluster1 = cluster[1], cluster2 = cluster[2]))
  }
}


simest.gm <- function(dat, rule, ni = 24) # rule = 0,1,2,3
{
  if (rule == 0) {
    dd <- dat |>
      group_by(fundid) |>
      mutate(inisize = first(lagsize)) |>
      ungroup()
  }
  if (rule == 1) {
    dd <- dat |>
      filter(survival1 == 1) |>
      group_by(fundid) |>
      filter(n() >= ni) |>
      mutate(inisize = first(lagsize), age = log(n())) |>
      ungroup()
  }
  if (rule == 2) {
    dd <- dat |>
      filter(survival2 == 1) |>
      group_by(fundid) |>
      filter(n() >= ni) |>
      mutate(inisize = first(lagsize), age = log(n())) |>
      ungroup()
  }
  if (rule == 3) {
    dd <- dat |>
      filter(survival3 == 1) |>
      group_by(fundid) |>
      filter(n() >= ni) |>
      mutate(inisize = first(lagsize), age = log(n())) |>
      ungroup()
  }

  dd_rd <- dd |>
    group_by(fundid) |>
    mutate(
      ret_RD = ret - rev(cumsum(rev(ret))) / c(n():1),
      logx_RD = lagsize - rev(cumsum(rev(lagsize))) / c(n():1)
    ) |>
    filter(!(ret_RD == 0 & logx_RD == 0)) |>
    ungroup()

  ## model fitting
  ols <- lm(ret ~ lagsize, data = dd) # ols
  fe <- feols(ret ~ lagsize | fundid, cluster = c("fundid"), data = dd, notes = FALSE) ## fe
  rdout <- estimation.rd(dat = dd_rd) # rd
  re <- lmer(ret ~ lagsize + (1 | fundid), data = dd) # re

  # browser()
  ## CW estimator
  sizefe <- feols(size ~ lagsize + I(ret + Bret) + ps | fundid, cluster = c("fundid"), data = dd, notes = FALSE)
  fe_fund <- fixef(sizefe)$fundid
  dd$fe_fund <- rep(fe_fund, table(dd$fundid))
  cw <- lmer(ret ~ lagsize + fe_fund + inisize + (1 | fundid), data = dd)
  cw2 <- lmer(ret ~ lagsize + fe_fund + (1 | fundid), data = dd)

  ## output
  best <- c(coef(ols)[2], coef(fe), rdout[1], fixef(re)[2], fixef(cw)[2], fixef(cw2)[2])

  test.ols <- coeftest.cluster(data = dd, fm = ols, cluster1 = c("fundid"))
  test.re <- robust.se.remodel(dat = dd, model = re, cluster = c("fundid"))
  test.cw <- robust.se.remodel(dat = dd, model = cw, cluster = c("fundid"))
  test.cw2 <- robust.se.remodel(dat = dd, model = cw2, cluster = c("fundid"))

  tv <- c(test.ols[2, 3], fe$coeftable[, 3], rdout[2], test.re[2, 3], test.cw[2, 3], test.cw2[2, 3])

  return(list(best = best, tv = tv))
}

# choose num of sim
nsim <- 1000

base_spec <- data.table(b = -0.0005, cor_a_stde = 0.15, cor_a_phi = 0.25, n_funds = 1000, n_periods = 120, note = "base")

b_alternative <- c(0, -.001, -.002) # c(0, -.001, -.002, -.003)
cor_a_stde_alternative <- c(0, 0.5)

n_periods_alternative <- c(60, 180, 240)
n_funds_alternative <- c(500, 2000)
cor_a_phi_alternative <- c(0, 0.125, 0.375, 0.5)

# b_alternative <- c(0)
# cor_a_stde_alternative <- c(0)

# we are going to vary one at a time
specs <- copy(base_spec)

for (this_b in b_alternative) {
  this_spec <- base_spec %>% mutate(b = this_b, note = "vary_b")
  specs <- rbind(specs, this_spec)
}

for (this_cor_a_stde in cor_a_stde_alternative) {
  this_spec <- base_spec %>% mutate(cor_a_stde = this_cor_a_stde, note = "vary_cor_a_stde")
  specs <- rbind(specs, this_spec)
}

for (this_n_periods in n_periods_alternative) {
  this_spec <- base_spec %>% mutate(n_periods = this_n_periods, note = "vary_n_periods")
  specs <- rbind(specs, this_spec)
}

for (this_n_funds in n_funds_alternative) {
  this_spec <- base_spec %>% mutate(n_funds = this_n_funds, note = "vary_n_funds")
  specs <- rbind(specs, this_spec)
}

for (this_cor_a_phi in cor_a_phi_alternative) {
  this_spec <- base_spec %>% mutate(cor_a_phi = this_cor_a_phi, note = "var_cor_a_phi")
  specs <- rbind(specs, this_spec)
}

specs[, spec_idx := .I]

# nb <- length(bseq) # number of b
nest <- 6 # number of estimators. Each fund exiting scenario has 6 estimators
ndropout <- 4 # number of fund exiting scenarios

## simulation PAUSE
to_dir <- paste0("../tmp/simulation/drs/nsim_", nsim, "/")
dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)

for (this_spec_idx in specs$spec_idx) {
  # this_spec_idx <- 1

  plan(multisession, workers = parallel::detectCores() - 2)

  b <- specs[spec_idx == this_spec_idx]$b
  cor_a_phi <- specs[spec_idx == this_spec_idx]$cor_a_phi
  cor_a_stde <- specs[spec_idx == this_spec_idx]$cor_a_stde
  n_funds <- specs[spec_idx == this_spec_idx]$n_funds
  n_periods <- specs[spec_idx == this_spec_idx]$n_periods

  # Wrap the simulation logic into a function

  sim_one_iteration <- function(i, b, cor_a_phi, cor_a_stde, n_funds, n_periods) {
    # 1. Generate data
    tmp <- gendata(
      gm = cor_a_phi, im = 0, b = b, N = n_funds, tt = n_periods,
      phi.mean = 0, phi.sd = 0.03, cor_a_stde = cor_a_stde
    )
    dd <- tmp$data

    # 2. Estimate DRS coefficients
    res0 <- simest.gm(dat = dd, rule = 0)
    res1 <- simest.gm(dat = dd, rule = 1)
    res2 <- simest.gm(dat = dd, rule = 2)
    res3 <- simest.gm(dat = dd, rule = 3)

    # 3. Return a list containing both result vectors
    return(list(
      best = c(res0$best, res1$best, res2$best, res3$best),
      btv  = c(res0$tv, res1$tv, res2$tv, res3$tv)
    ))
  }

  # Run in parallel
  # future.seed = TRUE is mandatory for valid random numbers in parallel
  tic(paste0("parallel simulation: spec_idx = ", this_spec_idx))

  # results <- sim_one_iteration(1, b, cor_a_phi, cor_a_stde, n_funds, n_periods)

  results <- future_lapply(1:nsim, sim_one_iteration,
    b = b, cor_a_phi = cor_a_phi, cor_a_stde = cor_a_stde, n_funds = n_funds, n_periods = n_periods,
    future.seed = TRUE,
    future.packages = c("data.table", "lme4")
  )
  toc()

  # Combine the results back into your matrices
  best <- do.call(rbind, lapply(results, `[[`, "best"))
  btv <- do.call(rbind, lapply(results, `[[`, "btv"))

  # Revert to sequential if done
  plan(sequential)

  bias <- (apply(best, 2, mean) - b) * 1000
  ss <- apply(best, 2, sd)
  rmse <- sqrt((bias / 1000)^2 + ss^2)
  frac.reject <- apply(btv, 2, FUN = function(x) sum(abs(x) > 1.96) / nsim)
  output <- data.table(
    b, cor_a_phi, cor_a_stde, n_funds, n_periods, nsim,
    method = rep(1:nest, ndropout),
    exit = rep(1:ndropout, nest), bias = bias, stderr = ss, rmse = rmse, rejectfrac = frac.reject
  )

  # parse method
  methods <- data.table(method = 1:nest, method_name = c("OLS", "FE", "RD", "RE", "CW1", "CW2"))
  output <- output[methods, on = .(method)]

  exits <- data.table(exit = 1:ndropout, exit_name = c("No fund exit", "Alpha-based fund exit", "Size-based fund exit", "Random fund exit"))
  output <- output[exits, on = .(exit)]

  saveRDS(output, paste0(to_dir, "spec_idx_", this_spec_idx, ".RDS"))
}

# # add names
# from_dir <- '../tmp/simulation/drs/nsim_1000/'
# files <- list.files(from_dir)
# for (this_file in files){
#   # this_file <- files[1]
#   print(this_file)
#   output <- readRDS(paste0(from_dir, this_file))
#   output[, exit := rep(1:4,6)]
#   exits <- data.table(exit = 1:4, exit_name = c("No fund exit", "Alpha-based fund exit", "Size-based fund exit", "Random fund exit"))
#   output <- output[exits, on = .(exit)]
#   saveRDS(output, paste0(from_dir, this_file))
# }
