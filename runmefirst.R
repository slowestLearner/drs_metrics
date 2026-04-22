rm(list = setdiff(ls(), "data_bk"))
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

# -- load functions
source("../../../data_M/cluster2.R")



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
