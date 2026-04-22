# --- calibrate survival probabilities out to 20 years
library(this.path)
setwd(this.path::this.dir())
source("../utility_functions/runmefirst.R")

data <- readRDS("../../../data/funds/min_sample.RDS")
data[, yyyy := floor(yyyymm / 100)]
data <- unique(data[, .(yyyy, fundid)])

max_y <- data[, max(yyyy)]
data <- data[, .(first_y = min(yyyy), last_y = max(yyyy)), fundid]
data[, max_possible_inferred_life_length := max_y - first_y]
data[, life_length := ifelse(last_y < max_y, last_y - first_y, 100)]

# cumulative death prob
p.cum_death_prob <- function(this_hor) {
  tmp <- data[max_possible_inferred_life_length >= this_hor]
  return(data.table(
    hor = this_hor, obs = nrow(tmp),
    cum_death_prob = tmp[, mean(life_length <= this_hor)]
  ))
}

out <- rbindlist(lapply(1:20, p.cum_death_prob))
to_dir <- "../tmp/simulation/inputs/"
dir.create(to_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(out, paste0(to_dir, "cum_death_prob.RDS"))

# plot(out[, .(hor, cum_death_prob)], type = "l")
