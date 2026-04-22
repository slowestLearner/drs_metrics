# --- Let me use Table 4 of Zhu (2018)
library(this.path)
setwd(this.path::this.dir())
source("../utility_functions/runmefirst.R")

# performance
ret <- c(0, 0, 0, -.03, -.2)
n <- length(ret)
ret_fwd_demeaned <- ret - rev(cumsum(rev(ret))) / c(n:1)
ret_fwd_demeaned <- ret_fwd_demeaned[-n]

data <- data.table(type = "raw return", idx = 1:n, ret = ret)
data <- rbind(data, data.table(type = "forward demeaned return", idx = 1:n, ret = c(ret_fwd_demeaned, NA)))


pp <- ggplot(data, aes(x = idx, y = ret, color = type)) +
  geom_line(lwd = 1.5) +
  geom_point(cex = 3) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(limits = factor(1:n)) +
  # 2. Use this logic for the legend
  scale_color_manual(
    values = c("blue", "red"), # Assign colors manually
    labels = parse(text = c("paste('Forward demeaned performance ', (bar(r)[t]))", "paste('Raw performance ', (r[t]))"))
  ) +
  theme(
    text = element_text(size = 15),
    legend.position = c(.3, .2),
    legend.title = element_blank()
  ) +
  labs(x = "Period", y = "Performance")

to_dir <- "../figs/bias/"
dir.create(to_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(paste0(to_dir, "illustration.png"), pp, device = "png", w = 6.5, h = 4.5, dpi = 600)
