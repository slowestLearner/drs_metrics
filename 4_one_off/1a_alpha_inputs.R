# --- Put together inputs (mu, Sigma) for subsequent portfolio optimization
library(this.path)
setwd(this.path::this.dir())
source("../runmefirst.R")

# get covariance and returns by sample
data <- readRDS("../../../data/factors/2x4x4_ff_portfolios.RDS")
vars <- setdiff(names(data), c("ret_vw", "ret_ew"))
data <- melt(data, id.vars = vars, variable.name = "ret_type", value.name = "ret") %>% mutate(ret_type = as.character(ret_type))

# mark the portfolios
tmp <- unique(data[, .(file, port)])[order(file)]
tmp[, port_id := .I]
tmp[, port_id := rank(port_id), file]
data <- data[tmp, on = .(file, port)]
data <- data[, .(yyyymm, file, port, port_id, ret_type, ret, me = num_firms * avg_size)]
rm(tmp)

# also get two subsamples
tmp <- unique(data[, .(yyyymm)])
tmp[, yyyy := floor(yyyymm / 100)]
tt <- unique(tmp[, .(yyyy)])[, bin_period := ntile(yyyy, 2)]
tt[, lab_period := paste0(min(yyyy), "-", max(yyyy)), bin_period]
tmp <- tmp[tt, on = .(yyyy)][, yyyy := NULL]
rm(tt)

data_alt <- copy(data)
data_alt <- data_alt[tmp, on = .(yyyymm)]

data[, bin_period := 0]
data[, lab_period := "full sample"]
stopifnot(dim(data) == dim(data_alt))
data <- rbindlist(list(data, data_alt))
rm(data_alt)

# compute means, covariances, and benchmark weights
p.get_one <- function(this_data) {
  # this_data <- data_list[[1]]

  # average return
  mu_data <- this_data[, .(mu = mean(ret)), .(bin_period, file, port_id, ret_type)]

  # covariance matrix
  ret_data <- dcast(this_data, yyyymm ~ port_id, value.var = "ret")[, yyyymm := NULL]
  cov_data <- cov(ret_data)
  cov_data <- data.table(port1_id = sort(unique(this_data$port_id)), as.data.table(cov_data))
  cov_data <- melt(cov_data, id.vars = "port1_id", variable.name = "port2_id", value.name = "cov")
  cov_data[, port2_id := as.integer(port2_id)]
  cov_data[, bin_period := this_data$bin_period[1]]
  cov_data[, file := this_data$file[1]]
  cov_data[, ret_type := this_data$ret_type[1]]

  # benchmark weights
  this_data[, w := me / sum(me), yyyymm]
  w_data <- this_data[, .(w = mean(w)), .(bin_period, file, port_id, ret_type)]
  stopifnot(abs(1 - w_data[, sum(w)]) < 1e-6)

  return(list(mu_data = mu_data, cov_data = cov_data, w_data = w_data))
}

data_list <- split(data, by = c("bin_period", "file", "ret_type"))
results <- lapply(data_list, p.get_one)

mu_data <- do.call(rbind, lapply(results, `[[`, "mu_data"))
cov_data <- do.call(rbind, lapply(results, `[[`, "cov_data"))
w_data <- do.call(rbind, lapply(results, `[[`, "w_data"))

# get period names
names <- unique(data[, .(bin_period, lab_period)])
mu_data <- mu_data[names, on = .(bin_period)]
cov_data <- cov_data[names, on = .(bin_period)]
w_data <- w_data[names, on = .(bin_period)]

# get portfolio ids
port_names <- unique(data[, .(file, bin_period, port_id, port)])
mu_data <- mu_data[port_names, on = .(file, bin_period, port_id)]
cov_data <- cov_data[port_names[, .(file, bin_period, port1_id = port_id, port1 = port)], on = .(file, bin_period, port1_id)]
cov_data <- cov_data[port_names[, .(file, bin_period, port2_id = port_id, port2 = port)], on = .(file, bin_period, port2_id)]
w_data <- w_data[port_names, on = .(file, bin_period, port_id)]

to_dir <- "../tmp/one_off/short_constrained_alpha/inputs/"
dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(mu_data, paste0(to_dir, "port_expected_ret.RDS"))
saveRDS(cov_data, paste0(to_dir, "port_covariance.RDS"))
saveRDS(w_data, paste0(to_dir, "port_benchmark_weights.RDS"))
