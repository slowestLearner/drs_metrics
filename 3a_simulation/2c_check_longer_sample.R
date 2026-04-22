library(this.path)
setwd(this.path::this.dir())
library(future.apply)
source("../utility_functions/runmefirst.R")

# check results after changing parameters. Just bias
from_dir <- "tmp/sim_data_240_periods/nsim_100_no_rho2_rho3"
data <- rbindlist(lapply(list.files(from_dir, full.names = T), fread))

# mark things
tmp <- data.table(exit = 1:4, exit_lab = c("non", "alpha-based", "size-based", "random"))
data <- data[tmp, on = .(exit)]

tmp <- data.table(method = 1:4, method_lab = c("FE", "RD", "RE", "WB"))
data <- data[tmp, on = .(method)]

ggplot(data[exit == 1], aes(x = cor_a_phi, y = bias, color = method_lab)) +
  geom_line(aes(color = method_lab), lwd = 3) +
  geom_point(aes(color = method_lab), cex = 5) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(text = element_text(size = 30), legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("no exit, no rho2 or rho3") +
  coord_cartesian(ylim = c(-.1, .05))
