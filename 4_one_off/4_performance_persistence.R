# --- Estimates yearly size-growth sensitivity. Used in the paper for a back-of-the-envelope calculation
library(this.path)
setwd(this.path::this.dir())
source("../utility_functions/runmefirst.R")


# summarize by year
data <- readRDS("../../data/funds/min_sample.RDS")
data[, yyyy := floor(yyyymm / 100)]
data <- data[, .(
    obs = .N,
    benchmark_adj_gret = sum(benchmark_adj_gret),
    gross_return = sum(gross_return),
    logtna = log(last(tna))
), .(yyyy, fundid)]
data <- data[obs == 12][, obs := NULL] %>% na.omit()

data <- merge(data, data[, .(yyyy = yyyy + 1, fundid, logtna_1 = logtna)], by = c("yyyy", "fundid"))

feols(I(logtna - logtna_1) ~ benchmark_adj_gret | yyyy, data = data, cluster = "fundid")
feols(I(logtna - logtna_1) ~ gross_return | yyyy, data = data, cluster = "fundid")

# # let's sort
# data[, bin := ntile(gross_return, 20), yyyy]

# out <- data[, .(
#     gross_return = mean(gross_return),
#     logtna_growth = mean(logtna - logtna_1)
# ), .(yyyy, bin)]

# out <- out[, .(gross_return = mean(gross_return), logtna_growth = mean(logtna_growth)), .(bin)]

# ggplot(out, aes(x = gross_return, y = logtna_growth)) +
#     geom_line(lwd = 1) +
#     geom_point() +
#     theme_classic() +
#     labs(x = "Gross return bin", y = "Log TNA growth") +
#     theme(text = element_text(size = 35))
