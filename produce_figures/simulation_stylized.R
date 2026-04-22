# -- this is simple, can always do this
library(this.path)
setwd(this.path::this.dir())
source("../utility_functions/runmefirst.R")
options(width = 150)

# read some inputs
files <- list.files("../3a_simulation/tmp/sim_data_1000/", full.names = T)
data <- rbindlist(lapply(files, fread))
data[, rejectfrac := NULL]
rm(files)

# parse a bit
tmp <- unique(data[, .(method)])[order(method)]
tmp[, method_lab := c("FE", "RD", "RE", "WB")]
data <- data[tmp, on = .(method)]

tmp <- data.table(exit = 1:4)
tmp[, exit_lab := c("no", "alpha", "size", "random")]
data <- data[tmp, on = .(exit)]

# change scale TODO: fix upstream!
data[, bias := bias / 10]
data[, stderr := stderr * 1e2]
data[, rmse := rmse * 1e2]
data_all <- copy(data)

# y-axis scales
# yy_bias <- range(data[, bias])
# yy_stderr <- c(0, max(data[, stderr]))
# yy_rmse <- c(0, max(data[, rmse]))
yy_bias <- c(-0.05, .02)
yy_stderr <- c(0, .008)
yy_rmse <- c(0, .05)

# plot
to_dir <- "../figs/simulation/stylized/baseline_params/"

for (this_exit in unique(data_all[, exit_lab])) {
    # this_exit <- 'no'
    data <- data_all[exit_lab == this_exit]

    pp <- ggplot(data, aes(x = cor_a_phi, y = bias, group = method_lab, color = method_lab, pch = method_lab, lty = method_lab)) +
        geom_line(lwd = 1) +
        geom_point() +
        geom_hline(yintercept = 0, lty = 2) +
        ylim(yy_bias) +
        theme_classic() +
        theme(text = element_text(size = 80), legend.position = "right", legend.title = element_blank()) +
        labs(
            x = expression(cor(a[i], phi[i])), y = "Bias"
            # y = expression("Bias" ~ (phantom() %*% 10^{3}))
        )
    ggsave(paste0(to_dir, this_exit, "_bias.png"), pp, device = "png", width = 4.5, height = 3, dpi = 600)

    pp <- ggplot(data, aes(x = cor_a_phi, y = stderr, group = method_lab, color = method_lab, pch = method_lab, lty = method_lab)) +
        geom_line(lwd = 1) +
        geom_point() +
        geom_hline(yintercept = 0, lty = 2) +
        ylim(yy_stderr) +
        theme_classic() +
        theme(text = element_text(size = 80), legend.position = c(1.8, .23), legend.title = element_blank()) +
        labs(
            x = expression(cor(a[i], phi[i])),
            y = "SD"
            # y = expression("SD" ~ (phantom() %*% 10^{3}))
        )
    ggsave(paste0(to_dir, this_exit, "_stderr.png"), pp, device = "png", width = 4, height = 3, dpi = 600)

    pp <- ggplot(data, aes(x = cor_a_phi, y = rmse, group = method_lab, color = method_lab, pch = method_lab, lty = method_lab)) +
        geom_line(lwd = 1) +
        geom_point() +
        geom_hline(yintercept = 0, lty = 2) +
        ylim(yy_rmse) +
        theme_classic() +
        theme(text = element_text(size = 80), legend.position = c(1.8, .23), legend.title = element_blank()) +
        labs(
            x = expression(cor(a[i], phi[i])),
            y = "RMSE"
            # y = expression("RMSE" ~ (phantom() %*% 10^{3}))
        )
    ggsave(paste0(to_dir, this_exit, "_rmse.png"), pp, device = "png", width = 4, height = 3, dpi = 600)
}
