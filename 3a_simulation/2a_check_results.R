library(this.path)
setwd(this.path::this.dir())
library(future.apply)
source("~/.runmefirst")

# check results after changing parameters. Just bias
from_dir <- "tmp/sim_data_1000/"
data <- rbindlist(lapply(list.files(from_dir, full.names = T), fread))
data <- data[method == 3, .(exit, spec_idx = 1, cor_a_phi, type = "old", bias)]

from_dir <- "tmp/sim_data_more_size_range/100_phi_mean_is_0"
tmp <- rbindlist(lapply(list.files(from_dir, full.names = T), fread))
tmp <- tmp[b == -0.0005 & method == 3, .(exit, spec_idx = 2, cor_a_phi, type = "phi.mean = 0, phi.sd = 0.03, max(q0) = 20,000", bias)]
data <- rbind(data, tmp)


from_dir <- "tmp/sim_data_more_size_range/100_phi_mean_0_phi_sd_0.01"
tmp <- rbindlist(lapply(list.files(from_dir, full.names = T), fread))
tmp <- tmp[b == -0.0005 & method == 3, .(exit, spec_idx = 3, cor_a_phi, type = "phi.mean = 0, phi.sd = 0.01, max(q0) = 20,000", bias)]
data <- rbind(data, tmp)

from_dir <- "tmp/sim_data_more_size_range/100_phi_mean_0_phi_sd_0.01_max_q0_100"
tmp <- rbindlist(lapply(list.files(from_dir, full.names = T), fread))
tmp <- tmp[b == -0.0005 & method == 3, .(exit, spec_idx = 4, cor_a_phi, type = "phi.mean = 0, phi.sd = 0.01, max(q0) = 100", bias)]
data <- rbind(data, tmp)
rm(tmp, from_dir)

out <- data[, .(bias = mean(bias)), by = .(cor_a_phi, spec_idx, type)]

ggplot(out, aes(x = cor_a_phi, y = bias, color = reorder(type, spec_idx), pch = reorder(type, spec_idx), lty = reorder(type, spec_idx))) +
  geom_line(lwd = 3) +
  geom_point(cex = 5) +
  theme_minimal() +
  theme(text = element_text(size = 30), legend.position = c(.3, .8), legend.title = element_blank(), plot.title = element_text(hjust = .5)) +
  ggtitle("RE bias")
