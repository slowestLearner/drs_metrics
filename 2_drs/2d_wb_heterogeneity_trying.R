# --- WB: heterogeneity. Just do baseline specification
library(this.path)
setwd(this.path::this.dir())
source("~/.runmefirst")
library(tidyverse)
library(lubridate)
library(lme4)
library(fixest)

# --- differences in approaches

# 1. I sorted by period
# 2. I've forced them to share the same intercept, which is not general
# 3.

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

# make it more interpretable
dat[, benchmark_adj_gret := 100 * benchmark_adj_gret]


# --- size model (the first one), get phi

size_model <- feols(logtna ~ lag_logtna + net_return + expanding_ret + logage + ner | fundid + yyyymm, cluster = c("fundid", "year"), dat)
int <- fixef(size_model)$fundid
dat$phi <- rep(int, as.numeric(table(dat$fundid)))
dat_all <- copy(dat)
rm(int, size_model)

# --- 1) variation with turnover

# - by cross-sectional, borderline sig

dat <- copy(dat_all)[!is.na(turnover)]
dat[, bin := ntile(turnover, 5), yyyymm]

mod <- lmer(benchmark_adj_gret ~ lag_logtna * as.factor(bin) + phi + si + (1 | fundid), data = dat)
obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
obj

# - by fund, again borderline

dat <- copy(dat_all)
tmp <- dat[, .(turnover = mean(turnover, na.rm = TRUE)), fundid][!is.na(turnover)]
tmp[, bin := ntile(turnover, 4)][, turnover := NULL]
dat <- dat[tmp, on = .(fundid), nomatch = NULL]
dat <- dat[order(fundid, yyyymm)]
rm(tmp)

mod <- lmer(benchmark_adj_gret ~ lag_logtna * as.factor(bin) + phi + si + (1 | fundid), data = dat)
obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
obj

# --- 2) variation with tracking error (SIG)

# - by fund, again borderline

dat <- copy(dat_all)
tmp <- dat[, .(tracking_error = sd(benchmark_adj_gret)), fundid]
tmp[, bin := ntile(tracking_error, 5)]
dat <- dat[tmp, on = .(fundid), nomatch = NULL][order(fundid, yyyymm)]
rm(tmp)

dat[, mean(tracking_error), bin][order(bin)]

mod <- lmer(benchmark_adj_gret ~ lag_logtna * as.factor(bin) + phi + si + (1 | fundid), data = dat)
obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
obj

# --- 3) variation with expense ratio (SIG)

# - in cross-section. bin 5 is a bit of weird
dat <- copy(dat_all)[!is.na(ner)]
dat[, bin := ntile(ner, 4), yyyymm]

dat[, mean(ner * 12), bin][order(bin)]

mod <- lmer(benchmark_adj_gret ~ lag_logtna * as.factor(bin) + phi + si + (1 | fundid), data = dat)
obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
obj


# - by fund, better

dat <- copy(dat_all)
tmp <- dat[!is.na(ner), .(ner = mean(ner)), fundid]
tmp[, bin := ntile(ner, 4)][, ner := NULL]
dat <- dat[tmp, on = .(fundid), nomatch = NULL][order(fundid, yyyymm)]
rm(tmp)

dat[, mean(ner), bin][order(bin)]

mod <- lmer(benchmark_adj_gret ~ lag_logtna * as.factor(bin) + phi + si + (1 | fundid), data = dat)
obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
obj

# --- what about smb loading?

dat <- copy(dat_all)

tmp <- copy(dat[, .(yyyymm, fundid, gross_return)])
tt <- fread("../../../data/factors/ff5_plus_mom_monthly_factors.csv")[, .(yyyymm, mktRf, smb)]
tmp <- tmp[tt, on = .(yyyymm), nomatch = NULL]

tt <- tmp[, .(obs = .N), fundid][obs >= 24]

p.get_one_fund <- function(this_fundid) {
    # this_fundid <- tt[1, fundid]
    mm <- lm(gross_return ~ mktRf + smb, data = tmp[fundid == this_fundid])
    return(data.table(fundid = this_fundid, smb_loading = mm$coefficients[3]))
}


tt <- rbindlist(lapply(tt$fundid, p.get_one_fund))
dat <- dat[tt, on = .(fundid), nomatch = NULL]
dat_bk <- copy(dat)
rm(tt, tmp, p.get_one_fund)

# -- by period, borderline
dat <- copy(dat_bk)
dat[, bin := ntile(smb_loading, 4), yyyymm]

dat[, mean(smb_loading), bin][order(bin)]

mod <- lmer(benchmark_adj_gret ~ lag_logtna * as.factor(bin) + phi + si + (1 | fundid), data = dat)
obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
obj

# -- by fund
dat <- copy(dat_bk)
tt <- dat[, .(smb_loading = last(smb_loading)), fundid][, bin := ntile(smb_loading, 4)][, smb_loading := NULL]
dat <- dat[tt, on = .(fundid), nomatch = NULL][order(fundid, yyyymm)]

dat[, mean(smb_loading), bin][order(bin)]

mod <- lmer(benchmark_adj_gret ~ lag_logtna * as.factor(bin) + phi + si + (1 | fundid), data = dat)
obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
obj

# --- liquidity of underlying funds? no not clear, really not clear

dat <- copy(dat_all)
tt <- unique(dat[, .(morningstar_category)])[order(morningstar_category)]
tt[, bin := c(1, 1, 1, 2, 2, 2, 3, 3, 3)]
dat <- dat[tt, on = .(morningstar_category), nomatch = NULL]
rm(tt)

mod <- lmer(benchmark_adj_gret ~ lag_logtna * as.factor(bin) + phi + si + (1 | fundid), data = dat)
obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
obj


# no, zero evidence!
mod <- lmer(benchmark_adj_gret ~ lag_logtna_bin1 + lag_logtna_bin2 + lag_logtna_bin3 + phi + si + (1 | fundid), data = dat)
r2_val <- suppressWarnings(performance::r2(mod))
obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))


# --- finally, what about size scales?

dat <- copy(dat_all)
dat[, bin := ntile(logtna, 4)]
dat[, .(min_tna = exp(min(logtna)), max_tna = exp(max(logtna))), bin][order(bin)]

mod <- lmer(benchmark_adj_gret ~ lag_logtna * as.factor(bin) + phi + si + (1 | fundid), data = dat)
obj <- robust.se.remodel(dat = dat, model = mod, cluster = c("fundid", "year"))
obj
