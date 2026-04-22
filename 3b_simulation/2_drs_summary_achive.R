# --- code to summarize DRS and plot
library(this.path)
setwd(this.path::this.dir())
source("../runmefirst.R")
library(ggpubr) # to place multiple ggplots on one page

# read results
from_dir <- "../tmp/simulation/drs/nsim_1000/"
files <- list.files(from_dir)
data_all <- rbindlist(lapply(paste0(from_dir, files), readRDS))
data_all[, method_name := paste0(method, "_", method_name)]
data_all[, exit_name := paste0(exit, "_", exit_name)]
rm(from_dir, files)

# canonical specs
base_spec <- data.table(b = -0.0005, cor_a_stde = 0.15, cor_a_phi = 0.25, n_funds = 1000, n_periods = 120, note = "base")

# Let's plot a few things
data <- copy(data_all[exit == 2 &
  # method_name == "CW1" &
  cor_a_phi == base_spec$cor_a_phi &
  n_funds == base_spec$n_funds & cor_a_stde == base_spec$cor_a_stde &
  n_periods == base_spec$n_periods])

dcast(data[, .(b, method_name, bias = round(bias, 4))], b ~ method_name, value.var = "bias")

# -- no correlation case
data <- copy(data_all[exit == 2 & n_funds == base_spec$n_funds & cor_a_stde == base_spec$cor_a_stde &
  n_periods == base_spec$n_periods])









from_dir <- "../tmp/simulation/drs/100_nsim/"

### put together all results
data <- data.frame()
method <- c("OLS", "FE", "RD", "RE", "CW1", "CW2")
exit <- c(0, 1, 2, 3)

for (i in 1:3) {
  # i <- 1
  # if (i == 1) b <- 0
  # if (i == 2) b <- -0.5
  # if (i == 3) b <- -1
  b <- bseq[i]

  for (j in 1:3) {
    name_c <- paste0("cor", j)
    name <- paste0(paste(paste0("Newoutput_b", i), name_c, sep = "_"), ".csv")

    tmp <- fread(paste0(from_dir, name))
    # tmp <- read_excel(name, sheet = 1)
    # tmp <- read_excel(name, sheet = 2)

    mcor <- mcorseq[j]
    # if (j == 1) mcor <- 0
    # if (j == 2) mcor <- 0.25
    # if (j == 3) mcor <- 0.5

    tmp$method <- rep(method, length(exit))
    tmp$exit <- rep(exit, each = length(method))
    tmp$b <- b
    tmp$mcor <- mcor
    # cross_df <- expand.grid(v2, v1)
    # estimator <- paste0(cross_df$Var2, cross_df$Var1)

    data <- rbind(data, tmp)
    # surv = rbind(surv, cbind(ss,  b, mcor))
  }
}

# names(data) =c("Estimator","exit","bias", "stderr", "rmse", "b", "mcor")
# names(surv) =c("survrate1", "survrate2", "survrate3", "b", "mcor")


data[, V1 := NULL]
setnames(data, "method", "Estimator")

to_dir <- "../figs/simulation/min/"
dir.create(to_dir, showWarnings = FALSE, recursive = TRUE)

# function to create a row of plots
plotfn <- function(plotdat, title, legend = "bottom") {
  # Plot Bias
  p1 <- ggplot(plotdat, aes(x = mcor, y = bias, group = Estimator, color = Estimator, pch = Estimator, lty = Estimator)) +
    geom_line() +
    geom_point() +
    xlab(expression(cor(a[i], phi[i]))) +
    ylab(expression("Bias" ~ (phantom() %*% 10^{
      3
    }))) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5)) +
    # xlab("Matching Correlation")  , lty=Estimator
    # ggtitle("Bias") +
    # ylim(biasrange) +
    theme_minimal()

  # Plot Standard Error
  p2 <- ggplot(plotdat, aes(x = mcor, y = stderr, group = Estimator, color = Estimator, pch = Estimator, lty = Estimator)) +
    geom_line() +
    geom_point() +
    xlab(expression(cor(a[i], phi[i]))) +
    ylab(expression("SD" ~ (phantom() %*% 10^{
      3
    }))) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5)) +
    # labs(title = "Standard Error") +
    theme_minimal()

  # Plot RMSE
  p3 <- ggplot(plotdat, aes(x = mcor, y = rmse, group = Estimator, color = Estimator, pch = Estimator, lty = Estimator)) +
    geom_line() +
    geom_point() +
    xlab(expression(cor(a[i], phi[i]))) +
    ylab(expression("RMSE" ~ (phantom() %*% 10^{
      3
    }))) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5)) +
    # labs(title = "RMSE") +
    theme_minimal()

  # Combine plots and add legend
  p <- ggarrange(p1, p2, p3, common.legend = TRUE, legend = legend, ncol = 3, nrow = 1)
  pp <- annotate_figure(p, top = text_grob(title, color = "black", face = "bold", size = 14))

  return(pp)
}

data <- plotdat <- data |> filter(Estimator != "CW1")
data$Estimator[data$Estimator == "CW2"] <- "WB"

for (i in 1:3) {
  bv <- bseq[i]
  # if (i == 1) bv <- 0
  # if (i == 2) bv <- -0.5
  # if (i == 3) bv <- -1

  for (j in 0:3) {
    if (j == 0) titlev <- "No fund exit"
    if (j == 1) titlev <- "Alpha-based fund exit"
    if (j == 2) titlev <- "Size-based fund exit"
    if (j == 3) titlev <- "Random fund exit"

    name1 <- paste0(paste(paste0("figure_b", i), paste0("e", j), sep = "_"), ".pdf")
    # name2 = paste0(paste(paste0("figure_b", i), paste0("e",j), sep="_"), "_1.pdf")

    # plotdat = data |> filter(b==bv & exit==j ) |> mutate(stderr = stderr * 10^3, rmse=rmse*10^3)
    plotdat <- data |>
      filter(b == bv & exit == j & Estimator != "OLS") |>
      mutate(stderr = stderr * 10^3, rmse = rmse * 10^3)

    flag <- ifelse(j < 3, "none", "bottom")
    pp <- plotfn(plotdat, title = titlev, legend = flag)

    pdf(paste0(to_dir, name1), width = 10, height = 4)
    print(pp)
    dev.off()
  }
}


# --- here are some one-off summaries
from_dir <- "../tmp/simulation/drs/100_nsim/"
files <- list.files(from_dir)
