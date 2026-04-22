robust.se.remodel <- function(dat, model, cluster = c("fundid", "year"), weights = FALSE) {
  # 1. Get the row names from the model - this is the ONLY reliable key
  # It tells us exactly which rows in 'dat' were used (and in what order)
  used_rows <- row.names(model.frame(model))
  dat_used <- dat[as.numeric(used_rows), ]

  # 2. Extract components with drop=FALSE to prevent vector conversion
  # check.names=FALSE prevents "model.matrix.model.....1."
  dd <- data.frame(
    y = model.frame(model)[, 1],
    model.matrix(model)[, -1, drop = FALSE],
    check.names = FALSE
  )

  # 3. Match Random Effects by ID
  # This replaces the broken 'rep(..., table())' logic
  random_intercepts <- ranef(model)$fundid
  dd$raneff <- random_intercepts[as.character(dat_used$fundid), 1]

  # 4. De-mean the dependent variable
  dd$y <- dd$y - dd$raneff

  # Identify regressors
  xnames <- setdiff(names(dd), c("y", "raneff"))
  fmla <- as.formula(paste("y ~ ", paste(xnames, collapse = " + ")))

  if (weights) {
    # Ensure weights are aligned to the used observations
    dd$w <- dat_used$w
    m1 <- lm(fmla, data = dd, weights = w)
  } else {
    m1 <- lm(fmla, data = dd)
  }

  # 5. Attach cluster variables to 'dd' from the ALIGNED data
  dd$fundid <- dat_used$fundid
  dd$year <- dat_used$year
  if ("w" %in% names(dat_used)) dd$w <- dat_used$w

  # Call the cluster function
  if (length(cluster) == 1) {
    return(coeftest.cluster.re(data = dd, fm = m1, cluster1 = cluster[1], weights = weights))
  } else {
    return(coeftest.cluster.re(data = dd, fm = m1, cluster1 = cluster[1], cluster2 = cluster[2], weights = weights))
  }
}

coeftest.cluster.re <- function(data, fm, cluster1 = NULL, cluster2 = NULL, ret = "test", weights = FALSE, K0 = 0) {
  library(sandwich)
  library(lmtest)

  data <- as.data.frame(data)

  if (is.null(cluster1)) {
    vcov_mat <- vcovHC(fm, type = "HC0")
    return(if (ret == "cov") vcov_mat else coeftest(fm, vcov_mat))
  }

  # Re-calculate estimating functions: (resid + raneff) * X
  # This uses the data passed in which is already aligned
  if (weights) {
    est.fun <- (residuals(fm) + data$raneff) * model.matrix(fm) * data$w
  } else {
    est.fun <- (residuals(fm) + data$raneff) * model.matrix(fm)
  }

  # Clean NAs (though there shouldn't be any if robust.se.remodel is used)
  inc.obs <- complete.cases(est.fun)
  est.fun <- est.fun[inc.obs, , drop = FALSE]

  NROW <- nrow(est.fun)
  N <- length(residuals(fm))
  K <- fm$rank + K0

  # Core Covariance Logic
  get_cov <- function(cl_var) {
    cl_var <- factor(cl_var[inc.obs], exclude = NULL)

    # Sum estimating functions by cluster
    u <- apply(est.fun, 2, function(x) tapply(x, cl_var, sum))
    meat <- crossprod(u) / N

    M <- length(levels(cl_var))
    dfc <- (M / (M - 1)) * ((N - 1) / (N - K))

    return(dfc * (NROW / N) * sandwich(fm, meat = meat))
  }

  # Calculate first cluster
  cov1 <- get_cov(data[[cluster1]])

  if (is.null(cluster2)) {
    return(if (ret == "cov") cov1 else coeftest(fm, cov1))
  } else {
    # Two-way clustering
    cov2 <- get_cov(data[[cluster2]])
    cluster12 <- paste(data[[cluster1]], data[[cluster2]], sep = "_")
    cov12 <- get_cov(cluster12)

    covMCL <- cov1 + cov2 - cov12
    return(if (ret == "cov") covMCL else coeftest(fm, covMCL))
  }
}
