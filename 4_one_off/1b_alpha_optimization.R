# --- I probably need to use the first half to estimate and the second half to look at OOS
library(this.path)
setwd(this.path::this.dir())
source("../utility_functions/runmefirst.R")
library(Matrix)
library(osqp)

# get inputs
mu_data <- readRDS("../tmp/one_off/short_constrained_alpha/inputs/port_expected_ret.RDS")
cov_data <- readRDS("../tmp/one_off/short_constrained_alpha/inputs/port_covariance.RDS")
w_data <- readRDS("../tmp/one_off/short_constrained_alpha/inputs/port_benchmark_weights.RDS")

in_sample <- "1963-1994"
out_of_sample <- "1995-2025"

# specifications
specs <- unique(mu_data[, .(file, ret_type)])
specs[, spec_idx := .I]
mu_data <- mu_data[specs, on = .(file, ret_type)]
cov_data <- cov_data[specs, on = .(file, ret_type)]
w_data <- w_data[specs, on = .(file, ret_type)]

# do port optimization
solve_qp_osqp <- function(cov_matrix, mean_vector, benchmark_weights, gamma, lw_shrinkage = 0.05, tolerance = 1e-10) {
    # can do shrinkage
    n <- nrow(cov_matrix)
    value_diagonal <- mean(diag(cov_matrix))
    value_off_diagonal <- (sum(cov_matrix) - sum(diag(cov_matrix))) / (n * (n - 1))
    CT <- (value_diagonal - value_off_diagonal) * diag(n) + value_off_diagonal * matrix(rep(1, n^2), n)
    C <- (1 - lw_shrinkage) * cov_matrix + lw_shrinkage * CT

    # 1. FORCE Sparse Matrix types (The 'dgCMatrix' format OSQP requires)
    P <- gamma * C
    q <- as.numeric(-mean_vector)

    n <- length(q)

    # 2. Build A (Constraint matrix) - must also be sparse
    # x >= l_bound
    A <- diag(n)
    l <- as.numeric(-benchmark_weights)
    u <- rep(Inf, n)

    solver_settings <- osqpSettings(
        eps_abs  = 1e-14,
        eps_rel  = 1e-14,
        max_iter = 100000,
        verbose  = FALSE
    )

    # Now the solver should accept the arguments
    solver <- osqp(P, q, A, l, u, solver_settings)
    res <- solver$Solve()

    # achieved results?
    w_active <- as.matrix(res$x)
    w <- w_active + benchmark_weights

    return(data.table(gamma,
        active_port_ret = (t(w_active) %*% mean_vector)[1],
        active_port_vol = sqrt((t(w_active) %*% cov_matrix %*% w_active)[1, 1]),
        port_id = 1:nrow(cov_matrix),
        w_b = benchmark_weights,
        w_a = as.numeric(w_active)
    ))
}

# wrap around and extract vol
p.get_vol_squared_distance <- function(cov_matrix, mean_vector, benchmark_weights, gamma, lw_shrinkage = 0.05) {
    out <- solve_qp_osqp(cov_matrix, mean_vector, benchmark_weights, gamma, lw_shrinkage)
    return((out$active_port_vol[1] - .04 / sqrt(12))^2)
}

# get outperformance in sample
p.get_one_spec_in_sample <- function(this_spec_idx) {
    # this_spec_idx <- 1
    mean_vector <- mu_data[spec_idx == this_spec_idx & lab_period == in_sample][order(port_id)][, mu / 100]
    cov_matrix <- dcast(cov_data[spec_idx == this_spec_idx & lab_period == in_sample][order(port1_id, port2_id)], port1_id ~ port2_id, value.var = "cov")[, -1] / 10000
    cov_matrix <- as.matrix(cov_matrix)
    benchmark_weights <- w_data[spec_idx == this_spec_idx & lab_period == in_sample][order(port_id)][, w]
    stopifnot(abs(1 - sum(benchmark_weights)) < 1e-6)

    # Let's first tune gamma to get target active vol around 5%
    opt <- optimize(function(x) {
        p.get_vol_squared_distance(cov_matrix, mean_vector, benchmark_weights, gamma = x, lw_shrinkage = 0.05)
    }, interval = c(0, 100))

    target_gamma <- opt$minimum
    out <- solve_qp_osqp(cov_matrix, mean_vector, benchmark_weights, gamma = target_gamma)

    # evaluate OOS
    mean_vector_oos <- mu_data[spec_idx == this_spec_idx & lab_period == out_of_sample][order(port_id)][, mu / 100]
    cov_matrix_oos <- dcast(cov_data[spec_idx == this_spec_idx & lab_period == out_of_sample][order(port1_id, port2_id)], port1_id ~ port2_id, value.var = "cov")[, -1] / 10000
    cov_matrix_oos <- as.matrix(cov_matrix_oos)

    out <- out[order(port_id)]
    w_a <- as.matrix(out[, w_a])

    out[, active_port_ret_oos := (t(w_a) %*% mean_vector_oos)[1]]
    out[, active_port_vol_oos := sqrt((t(w_a) %*% cov_matrix_oos %*% w_a)[1, 1])]

    # get port names, etc
    port_names <- mu_data[spec_idx == this_spec_idx & lab_period == out_of_sample][, .(port_id, spec_idx, lab_period, file, port, ret_type)]
    out <- out[port_names, on = .(port_id), nomatch = NULL]
    return(out)
}

out <- rbindlist(lapply(unique(mu_data$spec_idx), p.get_one_spec_in_sample))

to_dir <- "../tmp/one_off/short_constrained_alpha/outputs/"
dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(out, paste0(to_dir, "outperformance_is_and_oos.RDS"))
