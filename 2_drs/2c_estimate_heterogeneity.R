# --- estimate DRS heterogeneity using full-sample regressions
library(this.path)
setwd(this.path::this.dir())
source("../utility_functions/runmefirst.R")

########################## load data  ##########################################

# data: merge with fund size, tc
data <- readRDS("../tmp/raw_data/fund_net_ret_risk_adjusted.RDS")
tmp <- readRDS("../../../data/funds/active_equity_fund_monthly_data.RDS")[, .(yyyymm, wficn, exp_ratio, lag_logtna = log(mtna_1), turnover = turn_ratio)]
tmp[, idx := frank(yyyymm, ties.method = "dense")]
tmp <- tmp[tmp[, .(idx = idx - 1, wficn, logtna = lag_logtna)], on = .(idx, wficn), nomatch = NULL][, idx := NULL]
data <- data[tmp, on = .(yyyymm, wficn), nomatch = NULL]
data <- data[exp_ratio < .05] # otherwise a bit weird
rm(tmp)

# sometimes cluster by year instead of year-month (seems to make little difference)
data[, year := floor(yyyymm / 100)]

# infer fund age from sample start
tmp <- unique(data[, .(yyyymm)])
tmp[, date := as.Date(as.character(yyyymm * 100 + 28), "%Y%m%d")]
data <- data[tmp, on = .(yyyymm), nomatch = NULL]
data[, min_date := min(date), by = wficn]
data[, fund_age := as.numeric(date - min_date) / 365][, min_date := NULL]
data[, logage := log(fund_age + 1)]
rm(tmp)

# add exp_ratio back into performance measures. Also turn these variables into percent
vars <- names(data)
vars <- vars[grepl("alpha|ret", vars)]
for (this_var in vars) {
  data[, (this_var) := (get(this_var) + exp_ratio / 12) * 100]
}
rm(this_var)

# sort
data <- data[order(yyyymm, wficn)]

# mark if a fund is dead by the end of the sample
max_yyyymm <- max(data$yyyymm)
data <- data |>
  group_by(wficn) |>
  mutate(dead = max(yyyymm) < max_yyyymm, logage = log(fund_age + 1), si = first(lag_logtna)) %>%
  setDT()
rm(max_yyyymm)


# --- sort by a few characteristics into bins in the cross-section

# estimate SMB loading
tmp <- fread("../../../data/factors/ff5_plus_mom_monthly_factors.csv")[, .(yyyymm, smb)]
data <- data[tmp, on = .(yyyymm), nomatch = NULL]
rm(tmp)

fund_smb <- data[, .(obs = .N, b_smb = cov(ret / 100, smb) / var(smb)), wficn] %>%
  select(-obs) %>%
  setDT()
fund_smb[is.na(b_smb), b_smb := median(fund_smb[, b_smb], na.rm = T)]
fund_smb[, b_smb := DescTools::Winsorize(b_smb, quantile(b_smb, probs = c(.01, .99)))]
data <- data[fund_smb, on = .(wficn), nomatch = NULL] %>%
  mutate(-smb) %>%
  setDT()
rm(fund_smb)
data_bk <- copy(data)


########################## full-sample, one-way interaction ##########################################

data <- copy(data_bk)
data[, bin_turnover := ntile(turnover, 4), yyyymm]
data[, bin_exp_ratio := ntile(exp_ratio, 4), yyyymm]
data[, bin_b_smb := ntile(b_smb, 4), yyyymm]

# operate on one output variable
p.get_one_var <- function(this_var = "ret") {
  fe <- fixef(sizefe)$wficn
  data$phi <- rep(fe, as.numeric(table(data$wficn)))
  data$wficn <- as.factor(data$wficn)

  # formula <- paste0(this_var, " ~ lag_logtna * as.factor(bin) + phi + si + (1 | wficn)") %>% as.formula()
  # mod <- lmer(formula, data = data)

  formula <- paste0(this_var, " ~ lag_logtna * as.factor(bin) + phi + si") %>% as.formula()
  mod <- lme(formula, data = data, random = ~ 1 | wficn)
  obj <- robust.se.remodel(data = data, model = mod, cluster = c("wficn", "year"))

  out <- data.table(
    perf_var = this_var,
    var = rownames(obj), coef = obj[, 1], tstat = obj[, 3]
  )
  out <- out[grepl("lag_logtna", var)]
  print(out)
}

# --- choose specification

sizefe <- feols(logtna ~ lag_logtna + ret + logage + exp_ratio | wficn + yyyymm, cluster = c("wficn", "yyyymm"), data = data)

# data[, bin := bin_turnover]
data[, bin := bin_exp_ratio]

# takes around a min
tic()
out <- rbindlist(mclapply(vars, p.get_one_var, mc.cores = min(length(vars), detectCores() - 2)))
toc()

# compile results
tmp <- unique(out[, .(perf_var)]) %>% mutate(perf_var_lab = paste0(row_number(), "_", perf_var))
out <- out[tmp, on = .(perf_var)][, perf_var := NULL]
out[, coef_round := round(coef, 3)]
out[, tstat_round := round(tstat, 3)]
rm(tmp)

# report results
options(width = 200)
dcast(out, perf_var_lab ~ var, value.var = "coef_round")
dcast(out, perf_var_lab ~ var, value.var = "tstat_round")

# --- let's double sort results

data <- copy(data_bk)
data[, bin_turnover := ntile(exp_ratio, 3), yyyymm]
data[, bin_b_smb := ntile(b_smb, 3), yyyymm]
data[, bin := paste0("turnover_", bin_turnover, "-b_smb_", bin_b_smb)]

sizefe <- feols(logtna ~ lag_logtna + ret + logage + exp_ratio | wficn + yyyymm, cluster = c("wficn", "yyyymm"), data = data)

# takes 2 mins a bit long
tic()
out <- rbindlist(mclapply(vars, p.get_one_var, mc.cores = min(length(vars), detectCores() - 2)))
out_all <- copy(out)
toc()

# parse bins
out[, var := gsub("lag_logtna", "", var)]
out[, var := gsub(".as.factor.bin.", "", var)]
out[var == "", var := "baseline"]

tmp <- unique(out[, .(var)])
tmp[, bin_turnover := c(rep(1, 3), rep(2, 3), rep(3, 3))]
tmp[, bin_b_smb := rep(1:3, 3)]
out <- out[tmp, on = .(var), nomatch = NULL]
rm(tmp)

out[, coef_round := round(coef, 3)]
out[, tstat_round := round(tstat, 3)]

this_perf_var <- "alpha_CAPM"
this_perf_var <- "alpha_FF6"
dcast(out[perf_var == this_perf_var], bin_turnover ~ bin_b_smb, value.var = "coef_round")
dcast(out[perf_var == this_perf_var], bin_turnover ~ bin_b_smb, value.var = "tstat_round")
