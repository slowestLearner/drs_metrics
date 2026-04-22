# --- visualize results
library(this.path)
setwd(this.path::this.dir())
source("../runmefirst.R")

# full sample stuff
data_all <- readRDS("../tmp/one_off/short_constrained_alpha/outputs/outperformance_is_and_oos.RDS")

# annualize
data <- unique(data_all[, .(file, ret_type, lab_period,
    active_port_ret_is = 12 * active_port_ret,
    active_port_vol_is = sqrt(12) * active_port_vol,
    active_port_ret_oos = 12 * active_port_ret_oos,
    active_port_vol_oos = sqrt(12) * active_port_vol_oos
)])

data

# choose one to look at


# parse port meaning
port_names <- readRDS("../../../data/factors/2x4x4_ff_portfolios.RDS")[, .(file, port, port_me, port2, port3, char2, char3)] %>% unique()
data <- data[port_names, on = .(file, port), nomatch = NULL]
rm(port_names)
