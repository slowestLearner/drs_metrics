# -- Let's try to visualize WB estimates across a bunch of specifications
library(this.path)
setwd(this.path::this.dir())
source("~/.runmefirst")
options(width = 150)

# inputs
data <- readRDS("../tmp/drs/estimation/wb.RDS")[var %in% c("lag_logtna", "phi") & size_model_num == 1]
data[, c("nobs", "marginalR2", "size_model_num") := NULL]
data <- data[ret.name != "capm_adj_gret"]
data[, coef := 100 * coef]
data[, se := 100 * se]

# order them and rename
tmp <- unique(data[, .(ret.name)])[, ret.name_idx := .I]
tmp[, ret.name_lab := c("Benchmark-adj", "FF3 alpha", "Carhart alpha", "FF5 alpha", "FF5M alpha")]
data <- data[tmp, on = .(ret.name)]
rm(tmp)
data_all <- copy(data)

for (this_var in unique(data_all[, var])) {
    print(this_var)
    # this_var <- 'lag_logtna'
    # this_var <- 'phi'
    data <- copy(data_all[var == this_var])

    yy <- c(min(data[, coef - 2 * se]), max(data[, coef + 2 * se]))
    if (this_var == "phi") {
        yy <- c(0, yy[2])
    } else if (this_var == "lag_logtna") {
        yy <- c(yy[1], 0)
    }

    to_dir <- paste0("../figs/empirical/wb/varying_specifications/", this_var, "/")
    dir.create(to_dir, showWarnings = FALSE, recursive = TRUE)

    for (this_size_type in unique(data[, size_type])) {
        # this_size_type <- 'logtna'
        # this_size_type <- 'logtna_mktcap'
        subdata <- data[size_type == this_size_type]
        pp <- ggplot(subdata, aes(x = ret.name_idx, y = coef)) +
            geom_bar(stat = "identity") +
            geom_errorbar(aes(ymin = coef - 1.96 * se, ymax = coef + 1.96 * se), width = .5, lwd = .5) +
            theme_classic() +
            scale_x_discrete(
                limits = subdata[var == this_var, factor(ret.name_idx)],
                labels = subdata[var == this_var, ret.name_lab]
            ) +
            theme(text = element_text(size = 30), axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = element_blank(), y = "Estimate") +
            coord_cartesian(ylim = yy)
        ggsave(paste0(to_dir, this_size_type, ".png"), pp, device = "png", width = 3, height = 3, dpi = 300)
    }
}
