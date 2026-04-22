# --- This code moves some of the input data from J's computer to the dropbox folder. No need to run this. This is just an attempt to keep track of where data came from. When we produce a replication codebase, I will change the code to download from WRDS, etc.

# load libraries
library(this.path)
setwd(this.path::this.dir())
source("../utility_functions/runmefirst.R")
library(zoo)

# -------- market cap ---------
data <- readRDS("../../../../../../Desktop/J-Leaves/data/stockprices/raw/monthly/msf_cleaned/2024.RDS")
data <- data[shrcd %in% 10:12, .(total_mktcap = sum(prc * shrout / 1e3, na.rm = T)), .(yyyymm)]
saveRDS(data, "../../../data/macro/total_mkt_cap.RDS")

# ------- 2x4x4 Fama-French portfolios -------

data <- readRDS("../../../../../../Desktop/J-Leaves/data/portfolios/tripleSorts/2x4x4/combined.RDS")
saveRDS(data, "../../../data/factors/2x4x4_ff_portfolios.RDS")

# ------- active fund samples -------

#  monthly fund returns
data <- readRDS("../../../../../../Desktop/J-Leaves/data/institutions/crspMutual/processed/fund_flows/all_fund_flows/monthly_ret_by_wficn/202405.RDS")[, obs_crsp_fundno := NULL]

# filter out ETFs and non-active funds
tmp <- readRDS("../../../../../../Desktop/J-Leaves/data/institutions/crspMutual/raw/fund_indicative_filled/by_wficn/202405.RDS")
tmp <- tmp[et_flag == "" & obj == "ED" & index_fund_flag == ""]
tmp[, c("et_flag", "obj", "index_fund_flag") := NULL]

# must have most investment in equities, and not too levered
tmp <- tmp[(per_com >= 50) & (per_com < 120)]

# cannot have too much other stuff
tmp <- tmp[per_bond < 30 & per_corp < 30]

# finally, get rid of names with index in it
tmp[, name_lower := tolower(fund_name_largest_share_class)]
tmp <- tmp[!grepl("index", name_lower)]
tmp <- tmp[, .(yyyymm, wficn, fund_name_largest_share_class, exp_ratio, turn_ratio, keep = 1)]
data <- tmp[data, on = .(yyyymm, wficn)]
rm(tmp)

# can fill forward and backward
data <- data[order(wficn, yyyymm)]
vars <- setdiff(names(data), c("yyyymm", "wficn", "mtna_1", "ret"))
for (this_v in vars) {
    setnames(data, this_v, "xx")
    data <- data[order(wficn, -yyyymm)]
    data[, xx := na.locf(xx, na.rm = F), wficn]
    data <- data[order(wficn, yyyymm)]
    data[, xx := na.locf(xx, na.rm = F), wficn]
    setnames(data, "xx", this_v)
}
gc()

data <- data %>% na.omit()
data[, keep := NULL]
data <- data[mtna_1 >= 1]

to_dir <- "../../../data/funds/"
dir.create(to_dir, recursive = T, showWarnings = F)
saveRDS(data, paste0(to_dir, "active_equity_fund_monthly_data.RDS"))
