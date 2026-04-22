# --- later directly download from CRSP. For now, just move
library(this.path)
setwd(this.path::this.dir())
source("../runmefirst.R")

load('../../../Data_M/MFdata.RData')
data <- data.table(dat)
data[, yyyymm := 100*year + month]
data[, c('year','month','caldt_end') := NULL]
saveRDS(data, '../../../data/funds/min_sample.RDS')
