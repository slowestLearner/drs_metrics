# --- directly modified from Min's code, Analysis_DRS_forJ.R
library(this.path)
setwd(this.path::this.dir())
source("../utility_functions/runmefirst.R")
library(tidyverse)
library(lubridate)
library(lme4)
library(fixest)

# basic functions
source("../cluster2.R")

# load data
dat <- readRDS("../tmp/raw_data/reg_table_min.RDS")
dat <- dat[order(fundid, yyyymm)]

dat <- dat[yyyymm %in% 199501:201412]

# estimate using various risk-adjusted dynamics (later parallelize)
ret_names <- c("benchmark_adj_gret", "capm_adj_gret", "ff3_adj_gret", "carhart_adj_gret", "ff5_adj_gret", "ff6_adj_gret")

tic()
plan(multisession, workers = parallel::detectCores() - 2)
out <- rbindlist(future_lapply(ret_names, function(x) p.drs(dat, x), future.seed = 123, future.packages = c("data.table")))
plan(sequential)
toc()

out[ret.name == "benchmark_adj_gret"]

to_dir <- "../tmp/drs/estimation/"
dir.create(to_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(out, paste0(to_dir, "existing_methods.RDS"))
