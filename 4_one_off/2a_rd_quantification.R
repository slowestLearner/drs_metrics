# --- Uses tabulated Zhu (2018) size-effect estimates to plot linear/log specifications and numerically integrate total decline.
library(this.path)
setwd(this.path::this.dir())
source("../utility_functions/runmefirst.R")

# Table 4 from Zhu (2018)
data <- data.table(
  size = c(20, 41, 71, 112, 174, 270, 432, 689, 1266, 3828),
  linear_spec = c(33.5, 24.9, 21.3, 16.8, 26.5, 18.8, 26.7, 25.1, 21.0, 7.3),
  log_spec = c(21.0, 16.9, 14.5, 12.7, 11.6, 10.7, 10.1, 8.8, 8.3, 7.9)
)
data <- melt(data, id.vars = "size", variable.name = "type", value.name = "effect") %>% mutate(type = as.character(type))
data_all <- copy(data)

# ggplot(data, aes(x = size, y = 12 * effect / 100, color = type)) +
#   geom_line() +
#   geom_point() +
#   geom_hline(yintercept = 0) +
#   scale_x_continuous(trans = "log10") +
#   theme(text = element_text(size = 35))

# integrate effects
data <- copy(data_all[type == "log_spec"])
diff(range(data[, log(size)]))

# Calculate the total decline (integral) for each specification
integration_results <- data_all %>%
  group_by(type) %>%
  summarise(
    # x is log2(size)
    # The integral is the sum of trapezoid areas: (y1+y2)/2 * (x2-x1)
    total_decline_bp = {
      x <- log2(size)
      y <- effect
      n <- length(x)
      sum(((y[-n] + y[-1]) / 2) * diff(x))
    }
  )
