# --- Use existing methods to estimate DRS
# directly modified from Min's code, Analysis_DRS_forJ.R
library(this.path)
setwd(this.path::this.dir())
source("../utility_functions/runmefirst.R")

# NOTE TO MIN: the code will NOT work if I do not load this, which overwrites some utility functions in runmefirst.R
# I have yet to figure out one unified set of utility functions that work for everything
source("../utility_functions/cluster2_archive.R")

# load data
dat <- readRDS("../tmp/raw_data/reg_table_min.RDS")[order(fundid, yyyymm)]

# # optional: restrict to Zhu (2018) sample period
# dat <- dat[yyyymm %in% 199501:201412]

# estimate using various risk-adjusted dynamics (later parallelize)
ret_names <- c("benchmark_adj_gret", "capm_adj_gret", "ff3_adj_gret", "carhart_adj_gret", "ff5_adj_gret", "ff6_adj_gret")

tic()
plan(multisession, workers = parallel::detectCores() - 2)
out <- rbindlist(future_lapply(ret_names, function(x) p.drs(dat, x), future.seed = 123, future.packages = c("data.table")))
plan(sequential)
toc()

# sve
to_file <- "../tmp/drs/estimation/existing_methods.RDS"
dir.create(dirname(to_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, to_file)
