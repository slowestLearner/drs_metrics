gc()
library(arrow)
library(fixest)
library(ggplot2)
library(parallel)
library(DescTools)
library(tictoc)
library(haven)
options(scipen = 10)
library(showtext)
library(sandwich)
library(dplyr)
library(this.path)
library(future.apply)
library(data.table)
library(tidyverse)
library(lubridate)
library(lme4)
library(nlme)

# # -- load functions
# source("../cluster2.R")

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



# Zhu (JFE 2018)
estimation.rd <- function(dat, clusters = c("fundid"), weights = FALSE) #  ZHU's method
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
    if (length(clusters) == 1) {
        # b1 <- coeftest.cluster(data = dat, fm = rd, cluster1 = "fundid")
        b1 <- coeftest.cluster(data = dat, fm = rd, cluster1 = clusters[1])
    } else if (length(clusters) == 2) {
        b1 <- coeftest.cluster(data = dat, fm = rd, cluster1 = clusters[1], cluster2 = clusters[2])
    }
    # b1 = coeftest.cluster(data=dat, fm=rd, cluster1="SectorMonth", cluster2="FundId")      # double Clustered


    res <- c(b1[, c(1, 3, 4)])
    return(res)
}



robust.se.remodel <- function(data, model, cluster = c("fundid", "year"), weights = FALSE) {
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
    dd$raneff <- rep(random_intercepts[, 1], as.numeric(table(data$fundid)))

    # De-mean the response variable by the random effects
    dd$y <- dd$y - dd$raneff

    xnames <- setdiff(names(dd), c("y", "raneff"))
    fmla <- as.formula(paste0("y ~ ", paste(xnames, collapse = " + ")))

    # Fit the OLS model
    if (weights) {
        # Assumes 'w' exists in 'dat' or 'dd'
        m1 <- lm(fmla, data = dd, weights = data$w)
    } else {
        m1 <- lm(fmla, data = dd)
    }

    dd$fundid <- data$fundid
    dd$year <- data$year

    # Standard error clustering
    if (length(cluster) == 1) {
        return(coeftest.cluster.re(data = dd, fm = m1, cluster1 = cluster[1], weights = weights))
    }
    if (length(cluster) == 2) {
        return(coeftest.cluster.re(data = dd, fm = m1, cluster1 = cluster[1], cluster2 = cluster[2], weights = weights))
    }
}


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
