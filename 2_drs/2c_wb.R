# --- Wb: main effects
library(this.path)
setwd(this.path::this.dir())
source("~/.runmefirst")
library(tidyverse)
library(lubridate)
library(lme4)
library(fixest)

# basic functions
source("../cluster2.R")

# load data
dat <- readRDS("../tmp/raw_data/reg_table_min.RDS")[order(fundid, yyyymm)]

# take out some variable never used
dat[, c(
    "name", "firm_name", "inception_date", "tna", "benchmark_adj_nret", "ger", "fam_tna", "peer_adj_gret", "peer_adj_nret", "caldt_bgn",
    "expanding_capmb", "expanding_carb", "expanding_ff3b", "expanding_ff5b", "expanding_ff6b"
) := NULL]

# a few definitions of size
common_vars <- names(dat)
common_vars <- common_vars[!grepl("laglogtna", common_vars)]
dat <- rbindlist(list(
    cbind(dat[, .(size_type = "logtna", logtna, lag_logtna)], dat[, ..common_vars]),
    cbind(dat[, .(size_type = "logtna_cpi", logtna = logtna_cpi, lag_logtna = lag_logtna_cpi)], dat[, ..common_vars]),
    cbind(dat[, .(size_type = "logtna_mktcap", logtna = logtna_mktcap, lag_logtna = lag_logtna_mktcap)], dat[, ..common_vars])
))


# --- save size model estimates
size_models <- list()

# Define your dimensions
size_types <- unique(dat[, size_type])
model_nums <- 1:4

for (st in size_types) {
    # st <- size_types[1]
    tic(st)
    for (m in model_nums) {
        # m <- model_nums[1]
        # Unique key for the list
        key <- paste0(st, "_", m)

        # Define your formula dynamically or with an if-switch
        if (m == 1) {
            formula <- paste0(st, " ~ lag_", st, " + net_return + expanding_ret + logage + ner | fundid + yyyymm")
        } else if (m == 2) {
            formula <- paste0(st, " ~ lag_", st, " + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid + yyyymm")
        } else if (m == 3) {
            formula <- paste0(st, " ~ lag_", st, " + net_return + expanding_ret + logage + ner | fundid")
        } else if (m == 4) {
            formula <- paste0(st, " ~ lag_", st, " + net_return + expanding_ret + logage + ner + expanding_ret_sd + expanding_bret + expanding_bret_sd | fundid")
        }

        size_models[[key]] <- feols(as.formula(formula), cluster = c("fundid", "year"), data = dat[size_type == st])
    }
    toc()
}

# save these models, btw
out <- data.table() # coefs
out_fe <- data.table() # fixed effects

for (st in size_types) {
    print(st)
    # st <- size_types[1]
    for (m in model_nums) {
        # m <- model_nums[1]

        # Unique key for the list
        key <- paste0(st, "_", m)

        out <- rbind(out, data.table(
            size_type = st, spec_idx = m, var = names(size_models[[key]]$coefficients),
            coef = size_models[[key]]$coefficients, se = size_models[[key]]$se,
            obs = size_models[[key]]$nobs, ar2 = r2(size_models[[key]])["ar2"], war2 = r2(size_models[[key]])["war2"]
        ))

        fe_all <- fixef(size_models[[key]])
        phi_vector <- fe_all$fundid

        out_fe <- rbind(out_fe, data.table(
            size_type = st, spec_idx = m,
            fundid = names(phi_vector),
            phi = as.numeric(phi_vector)
        ))
    }
}

to_dir <- "../tmp/drs/estimation/"
dir.create(to_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(out, paste0(to_dir, "size_dynamics.RDS"))
saveRDS(out_fe, paste0(to_dir, "size_dynamics_fe.RDS"))

## ----- WB estimates
# size_type <- 'logtna'
# size_model_num <- 1
p.wb <- function(dat = dat, size_models, size_type = "logtna", size_model_num, ret.name = "benchmark_adj_gret") {
    key <- paste0(size_type, "_", size_model_num)

    femod <- size_models[[key]]
    int <- fixef(femod)$fundid
    dat$phi <- rep(int, as.numeric(table(dat$fundid)))
    setnames(dat, ret.name, "lhs")

    mod <- lmer(lhs ~ lag_logtna + phi + si + (1 | fundid), data = dat)
    r2_val <- suppressWarnings(performance::r2(mod))
    obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))

    out <- data.table(
        size_type = size_type, size_model_num = size_model_num, var = rownames(obj),
        ret.name = ret.name,
        coef = obj[, 1], se = obj[, 2], nobs = nobs(mod), marginalR2 = unlist(r2_val[2])
    )

    return(out)
}

ret_types <- c("benchmark_adj_gret", "capm_adj_gret", "ff3_adj_gret", "carhart_adj_gret", "ff5_adj_gret", "ff6_adj_gret")

options(future.globals.maxSize = 2000 * 1024^2) # Set to 3GB

# this takes a while. each combination takes like 20-30 secs. there are 3 * 6 = 18 instances, so will take like 9 mins
out <- data.table()
for (this_size_type in size_types) {
    # this_size_type <- size_types[1]
    current_type_models <- size_models[paste0(this_size_type, "_", 1:4)]
    dat_type <- dat[size_type == this_size_type]

    for (ret_type in ret_types) {
        # ret_type <- ret_types[1]
        tic(paste0(this_size_type, "_", ret_type))
        plan(multisession, workers = parallel::detectCores() - 2)
        out <- rbind(out, rbindlist(future_lapply(1:4, function(i) p.wb(dat_type, current_type_models, this_size_type, i, ret_type), future.seed = 123, future.packages = c("data.table"))))
        plan(sequential)
        toc()
    }
}
out_all <- copy(out)
saveRDS(out_all, "../tmp/drs/estimation/wb.RDS")
gc()


# --- check variation of a_i and phi_i

this_size_type <- "logtna"
sub_dat <- copy(dat[size_type == this_size_type])

# estimate phi and coefs
i <- 1
femod <- size_models[[paste0(this_size_type, "_", i)]]
int <- fixef(femod)$fundid
sub_dat$phi <- rep(int, as.numeric(table(sub_dat$fundid)))
mod <- lmer(benchmark_adj_gret ~ lag_logtna + phi + si + (1 | fundid), data = sub_dat)
ss <- summary(mod)

sub_dat[, ret_no_size := benchmark_adj_gret - ss$coefficients[2, 1] * lag_logtna]
# sub_dat[, ret_no_size := benchmark_adj_gret + 0.0026 * lag_logtna] # alternative, zhu (2018) estimate
tmp <- sub_dat[, .(yyyymm, morningstar_category, fundid, ret_no_size, phi)]

# get rid of benchmark-based covariantion
tmp[, ret_no_size_demeaned := ret_no_size - mean(ret_no_size), yyyymm]
out <- tmp[, .(
    obs = .N, a = mean(ret_no_size_demeaned), se_a = sd(ret_no_size_demeaned) / sqrt(.N),
    phi = mean(phi), phi = last(phi)
), fundid]

# get summaries
p.summarize <- function(min_obs) {
    return(out[obs >= min_obs, .(
        min_obs = min_obs, num_funds = .N,
        sd_a_naive = sd(a),
        se_a = mean(se_a),
        sd_a_corrected = sqrt(sd(a)^2 - mean(se_a)^2),
        sd_phi = sd(phi), cor_a_phi = cor(a, phi)
    )])
}

options(width = 160)
rbindlist(lapply(c(0, 24, 36, 60, 120, 180, 240, 300), p.summarize))

# what about the expected gross return at a point in time?
sub_dat[, a := mean(ret_no_size), fundid]
sub_dat[, expected_gross_return := a - ss$coefficients[2, 1] * lag_logtna]
sub_dat[, sd(expected_gross_return), yyyymm][, mean(V1)]

sub_dat[, cor(a, lag_logtna), yyyymm][, mean(V1)]
sub_dat[, cor(a, lag_logtna, method = "spearman"), yyyymm][, mean(V1)]
