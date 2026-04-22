# R function for computing two-way cluster-robust standard errors.
# The code below was adapted by Ian Gow on 2011-05-16 using code supplied
# via Mitchell Petersen's website by Mahmood Arai, 2008-01-21. Modified
# on 2014-04-18 to return White (1980) standard errors if no cluster
# variable is provided and to add links to test code.
#
# Apart from a little cleanup of the code, the main difference between this
# and the earlier code is in the handling of missing values. Look at the file
# cluster.test.R to see example usage. Note that care should be taken to
# do subsetting outside of the call to lm or glm, as it is difficult to recon-
# struct subsetting of this kind from the fitted model. However, the code
# does handle transformations of variables in the model (e.g., logs). Please
# report any bugs, suggestions, or errors to iandgow@gmail.com.
#
# The output has  been tested fairly extensively against output of
# Mitchell Petersen'scluster2.ado commmand (hence implicitly against
# the Matlab and SAS code posted elsewhere here), but I have not tested
# the code against code for non-linear models, such as logit2.ado.
# For testing code, see here: http://www.iangow.me/~igow/code/cluster_test.Rnw
# For output, see here:http://www.iangow.me/~igow/code/cluster_test.pdf

# See: Thompson (2006), Cameron, Gelbach and Miller (2006) and Petersen (2010).
# and Gow, Ormazabal, and Taylor (2010) for more discussion of this code
# and two-way cluster-robust standard errors.

# The arguments of the function are data, fitted model, cluster1 and cluster2
# You need to install packages `sandwich' by Thomas Lumley and Achim Zeileis and
# `lmtest' by Torsten Hothorn, Achim Zeileis, Giovanni Millo and David Mitchell.
# (For example, type install.packages("sandwich") on the R console.)

coeftest.cluster <- function(data, fm, cluster1 = NULL, cluster2 = NULL, ret = "test", K0 = 0) {
  library(sandwich)
  library(lmtest)

  data <- as.data.frame(data)

  # Return White (1980) standard errors if no cluster
  # variable is provided
  if (is.null(cluster1)) {
    if (ret == "cov") {
      return(vcovHC(fm, type = "HC0"))
    } else {
      return(coeftest(fm, vcov = vcovHC(fm, type = "HC0")))
    }
  }

  # Calculation shared by covariance estimates
  # estimating function: (y - x'beta) * x
  est.fun <- estfun(fm)

  # Need to identify observations used in the regression (i.e.,
  # non-missing) values, as the cluster vectors come from the full
  # data set and may not be in the regression model.
  inc.obs <- !is.na(est.fun[, 1])
  est.fun <- est.fun[inc.obs, ]

  # print(length(est.fun))

  # Shared data for degrees-of-freedom corrections
  NROW <- NROW(est.fun)
  N <- dim(fm$model)[1]
  K <- fm$rank + K0

  # Calculate the sandwich covariance estimate
  cov <- function(cluster) {
    cluster <- factor(cluster, exclude = NULL)

    # Calculate the "meat" of the sandwich estimators
    if (length(est.fun) == length(cluster)) est.fun <- matrix(est.fun, ncol = 1)
    # print(length(est.fun))
    # print(length(cluster))

    u <- apply(est.fun, 2, function(x) tapply(x, cluster, sum))
    meat <- crossprod(u) / N

    # Calculations for degrees-of-freedom corrections, followed
    # by calculation of the variance-covariance estimate.
    # NOTE: NROW/N is a kluge to address the fact that sandwich
    # uses the wrong number of rows (includes rows omitted from
    # the regression).
    M <- length(levels(cluster))
    dfc <- M / (M - 1) * (N - 1) / (N - K)

    # print (sandwich(fm, meat=meat))
    return(dfc * NROW / N * sandwich(fm, meat = meat))
  }

  # Calculate the covariance matrix estimate for the first cluster.
  cluster1 <- data[inc.obs, cluster1]
  # browser()
  cov1 <- cov(cluster = cluster1)
  # print(cov1)

  if (is.null(cluster2)) {
    # If only one cluster supplied, return single cluster
    # results
    if (ret == "cov") {
      return(cov1)
    } else {
      return(coeftest(fm, cov1))
    }
  } else {
    # Otherwise do the calculations for the second cluster
    # and the "intersection" cluster.
    cluster2 <- data[inc.obs, cluster2]
    cluster12 <- paste(cluster1, cluster2, sep = "")

    # Calculate the covariance matrices for cluster2, the "intersection"
    # cluster, then then put all the pieces together.
    cov2 <- cov(cluster = cluster2)
    cov12 <- cov(cluster = cluster12)
    covMCL <- (cov1 + cov2 - cov12)

    # Return the output of coeftest using two-way cluster-robust
    # standard errors.
    # print(ret)
    if (ret == "cov") {
      return(covMCL)
    } else {
      return(coeftest(fm, covMCL))
    }
  }
}

# Following based on suggestion from
# https://stat.ethz.ch/pipermail/r-help/2011-January/264777.html
# provided by Achim Zeileis.
summary.cluster <- function(obj, data, cluster1, cluster2 = NULL, alpha = 0.05) {
  require(memisc)

  # Get original summary
  s <- getSummary(obj, alpha = alpha)

  ## replace Wald tests of coefficients
  s$coef[, 1:4] <- coeftest.cluster(data, obj, cluster1, cluster2)

  ## replace confidence intervals
  crit <- qt(alpha / 2, obj$df.residual)
  s$coef[, 5] <- s$coef[, 1] + crit * s$coef[, 2]
  s$coef[, 6] <- s$coef[, 1] - crit * s$coef[, 2]

  # Note that some components of s$sumsstat will be inconsistent with
  # the clustered calculations

  return(s)
}



coeftest.cluster.re <- function(data, fm, cluster1 = NULL, cluster2 = NULL, ret = "test", weights = FALSE, K0 = 0) {
  library(sandwich)
  library(lmtest)

  data <- as.data.frame(data)

  # Return White (1980) standard errors if no cluster
  # variable is provided
  if (is.null(cluster1)) {
    if (ret == "cov") {
      return(vcovHC(fm, type = "HC0"))
    } else {
      return(coeftest(fm, vcov = vcovHC(fm, type = "HC0")))
    }
  }

  # Calculation shared by covariance estimates
  # estimating function: (y - x'beta) * x
  # est.fun <- estfun(fm)
  if (weights) {
    est.fun <- (residuals(fm) + data$raneff) * model.matrix(fm) * data$w
  } else {
    est.fun <- (residuals(fm) + data$raneff) * model.matrix(fm)
  }
  # Need to identify observations used in the regression (i.e.,
  # non-missing) values, as the cluster vectors come from the full
  # data set and may not be in the regression model.
  inc.obs <- !is.na(est.fun[, 1])
  est.fun <- est.fun[inc.obs, ]

  # print(length(est.fun))

  # Shared data for degrees-of-freedom corrections
  NROW <- NROW(est.fun)
  N <- dim(fm$model)[1]
  K <- fm$rank + K0

  # Calculate the sandwich covariance estimate
  cov <- function(cluster) {
    cluster <- factor(cluster, exclude = NULL)

    # Calculate the "meat" of the sandwich estimators
    if (length(est.fun) == length(cluster)) est.fun <- matrix(est.fun, ncol = 1)
    # print(length(est.fun))
    # print(length(cluster))

    u <- apply(est.fun, 2, function(x) tapply(x, cluster, sum))
    meat <- crossprod(u) / N

    # Calculations for degrees-of-freedom corrections, followed
    # by calculation of the variance-covariance estimate.
    # NOTE: NROW/N is a kluge to address the fact that sandwich
    # uses the wrong number of rows (includes rows omitted from
    # the regression).
    M <- length(levels(cluster))
    dfc <- M / (M - 1) * (N - 1) / (N - K)

    # print (sandwich(fm, meat=meat))
    return(dfc * NROW / N * sandwich(fm, meat = meat))
  }

  # Calculate the covariance matrix estimate for the first cluster.
  cluster1 <- data[inc.obs, cluster1]
  # browser()
  cov1 <- cov(cluster = cluster1)
  # print(cov1)

  if (is.null(cluster2)) {
    # If only one cluster supplied, return single cluster
    # results
    if (ret == "cov") {
      return(cov1)
    } else {
      return(coeftest(fm, cov1))
    }
  } else {
    # Otherwise do the calculations for the second cluster
    # and the "intersection" cluster.
    cluster2 <- data[inc.obs, cluster2]
    cluster12 <- paste(cluster1, cluster2, sep = "")

    # Calculate the covariance matrices for cluster2, the "intersection"
    # cluster, then then put all the pieces together.
    cov2 <- cov(cluster = cluster2)
    cov12 <- cov(cluster = cluster12)
    covMCL <- (cov1 + cov2 - cov12)

    # Return the output of coeftest using two-way cluster-robust
    # standard errors.
    # print(ret)
    if (ret == "cov") {
      return(covMCL)
    } else {
      return(coeftest(fm, covMCL))
    }
  }
}


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
