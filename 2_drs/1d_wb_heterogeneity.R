# --- WB: heterogeneity. Sort by fund-level statistics.
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
dat <- readRDS("../tmp/raw_data/reg_table_min.RDS")[order(fundid, yyyymm)]

# take out some variable never used
dat[, c(
    "name", "firm_name", "inception_date", "tna", "benchmark_adj_nret", "ger", "fam_tna", "peer_adj_gret", "peer_adj_nret", "caldt_bgn",
    "expanding_capmb", "expanding_carb", "expanding_ff3b", "expanding_ff5b", "expanding_ff6b",
    "logtna_cpi", "logtna_mktcap",
    "lag_logtna_cpi", "lag_logtna_mktcap"
) := NULL]

# size model (the first one), get phi
size_model <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner | fundid + yyyymm, cluster = c("fundid", "year"), dat)
int <- fixef(size_model)$fundid
dat$phi <- rep(int, as.numeric(table(dat$fundid)))
rm(int, size_model)

# use SMB loadings to represent style
tmp <- copy(dat[, .(yyyymm, fundid, gross_return)])
tt <- fread("../../../data/factors/ff5_plus_mom_monthly_factors.csv")[, .(yyyymm, mktRf, smb)]
tmp <- tmp[tt, on = .(yyyymm), nomatch = NULL]
tt <- tmp[, .(obs = .N), fundid]

p.get_one_fund <- function(this_fundid) {
    # this_fundid <- tt[1, fundid]
    mm <- lm(gross_return ~ mktRf + smb, data = tmp[fundid == this_fundid])
    return(data.table(fundid = this_fundid, smb_loading = mm$coefficients[3]))
}

tt <- rbindlist(lapply(tt$fundid, p.get_one_fund))
dat <- merge(dat, tt, by = "fundid", all.x = T)
dat_all <- copy(dat)
rm(tt)

# --- summarize fund-level estimates

fund_summaries <- dat[, .(tracking_error = sd(benchmark_adj_gret), ner = mean(ner, na.rm = TRUE), smb_loading = mean(smb_loading, na.rm = TRUE)), fundid]

# vary a few specifications
ret_types <- c("benchmark_adj_gret", "ff3_adj_gret", "carhart_adj_gret", "ff5_adj_gret", "ff6_adj_gret")
sort_types <- c("tracking_error", "ner", "smb_loading")

out <- data.table()
out_sort_var <- data.table()
for (this_ret_type in ret_types) {
    # this_ret_type <- ret_types[1]
    print(this_ret_type)

    for (this_sort_type in sort_types) {
        # this_sort_type <- sort_types[1]
        dat <- copy(dat_all)

        sort_dat <- fund_summaries[, c("fundid", this_sort_type), with = F]
        setnames(sort_dat, this_sort_type, "sort_value")
        sort_dat <- sort_dat[!is.na(sort_value)]
        sort_dat[, bin := ntile(sort_value, 4)]

        if (this_sort_type == sort_types[1]) {
            out_sort_var <- rbind(out_sort_var, sort_dat[, .(sort_type = this_sort_type, sort_value = mean(sort_value)), bin])
        }


        dat <- merge(dat, sort_dat, by = "fundid")
        mod <- lmer(get(this_ret_type) ~ lag_logtna * as.factor(bin) + phi + si + (1 | fundid), data = dat)
        r2_val <- suppressWarnings(performance::r2(mod))
        obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
        out <- rbind(out, data.table(
            ret_type = this_ret_type, sort_type = this_sort_type,
            var = rownames(obj), coef = obj[, 1], se = obj[, 2], nobs = nobs(mod), marginalR2 = r2_val[[2]]
        ))
    }
}

to_dir <- "../tmp/drs/estimation/wb_hetero/"
dir.create(to_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(out, paste0(to_dir, "reg_results.RDS"))
saveRDS(out_sort_var, paste0(to_dir, "sort_var.RDS"))
