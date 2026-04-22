library(this.path)
setwd(this.path::this.dir())
library(future.apply)
source("../utility_functions/runmefirst.R")

# check results after changing parameters. Just bias
from_dir <- "tmp/1d/240_periods_with_new_WB/phi_mean_0.12_phi_sd_0.03/nsim_30_q0_max_100/"
data <- rbindlist(lapply(list.files(from_dir, full.names = T), fread))

# mark things
tmp <- data.table(exit = 1:4, exit_lab = c("none", "alpha-based", "size-based", "random"))
data <- data[tmp, on = .(exit)]

tmp <- data.table(method = 1:5, method_lab = c("FE", "RD", "RE", "WB", "new WB"))
data <- data[tmp, on = .(method), nomatch = NULL]

yy <- c(-1.3, .3)

this_exit <- 4
ggplot(data[exit == this_exit], aes(x = cor_a_phi, y = bias, color = reorder(method_lab, method))) +
  geom_line(aes(color = reorder(method_lab, method)), lwd = 1) +
  geom_point(aes(color = reorder(method_lab, method)), cex = 2) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(text = element_text(size = 11), legend.position = "bottom", legend.title = element_blank()) +
  coord_cartesian(ylim = yy) + 
  ggtitle(paste0('exit=', data[exit == this_exit][1, exit_lab]))

