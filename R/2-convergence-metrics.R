# Script info -------------------------------------------------------------

## File: 2-performance-metrics.R
## Owner: Matías Castillo-Aguilar
## This script performs:
##   1. RRi curve simulation using literature-based parameters from the Castillo-Aguilar et al. (2025) RRi-vs-time model.
##   2. Mapping of true parameters from original model to the reparameterized version of the model.
##   3. Narrow, moderate and wide Prior-specification from true parameters.
##   4. Fitting of prior-only Bayesian model to avoid future redundant model recompiling.
##   5. Model fitting for each simulated RRi curve using precompiled model.

# System info -------------------------------------------------------------

# sessionInfo()
#> R version 4.4.2 (2024-10-31 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 11 x64 (build 26100)
#>
#> Matrix products: default
#>
#>
#> locale:
#> [1] LC_COLLATE=Spanish_Chile.utf8  LC_CTYPE=Spanish_Chile.utf8    LC_MONETARY=Spanish_Chile.utf8
#> [4] LC_NUMERIC=C                   LC_TIME=Spanish_Chile.utf8
#>
#> time zone: America/Punta_Arenas
#> tzcode source: internal
#>
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base
#>
#> loaded via a namespace (and not attached):
#>  [1] bridgesampling_1.1-2 tensorA_0.36.2.1     future_1.49.0        generics_0.1.4       stringi_1.8.7
#>  [6] lattice_0.22-7       listenv_0.9.1        digest_0.6.37        magrittr_2.0.3       grid_4.4.2
#> [11] estimability_1.5.1   RColorBrewer_1.1-3   mvtnorm_1.3-3        plyr_1.8.9           Matrix_1.7-3
#> [16] pkgbuild_1.4.8       backports_1.5.0      gridExtra_2.3        Brobdingnag_1.2-9    purrr_1.0.4
#> [21] brms_2.22.0          QuickJSR_1.7.0       scales_1.4.0         codetools_0.2-20     abind_1.4-8
#> [26] cli_3.6.5            rlang_1.1.6          parallelly_1.44.0    future.apply_1.11.3  StanHeaders_2.32.10
#> [31] tools_4.4.2          rstan_2.32.7         inline_0.3.21        parallel_4.4.2       reshape2_1.4.4
#> [36] rstantools_2.4.0     checkmate_2.3.2      coda_0.19-4.1        dplyr_1.1.4          ggplot2_3.5.2
#> [41] globals_0.18.0       vctrs_0.6.5          posterior_1.6.1      R6_2.6.1             matrixStats_1.5.0
#> [46] stats4_4.4.2         lifecycle_1.0.4      emmeans_1.11.1       stringr_1.5.1        pkgconfig_2.0.3
#> [51] RcppParallel_5.1.10  pillar_1.10.2        gtable_0.3.6         loo_2.8.0            glue_1.8.0
#> [56] Rcpp_1.0.14          tibble_3.2.1         tidyselect_1.2.1     rstudioapi_0.17.1    farver_2.1.2
#> [61] xtable_1.8-4         bayesplot_1.12.0     nlme_3.1-168         compiler_4.4.2       distributional_0.5.0

# 0. Prepare workspace -------------------------------------------------------

## Load libraries
library(data.table)
library(brms)
library(ggplot2)

## Load custom functions
source("R/_functions.R")

## Identify model files for future loading
newFiles <- list.files(path = "models/new/sim", pattern = "\\.RDS", full.names = TRUE)
oldFiles <- list.files(path = "models/old/sim", pattern = "\\.RDS", full.names = TRUE)

# Estimating convergence measures -----------------------------------------

## Avoid computing if file already exists
if (!file.exists("output/rawConvergence.RDS")) {
  ## For new model
  newConvergence <- lapply(newFiles, function(i) {
    model <- readRDS(i)$fit
    list(
      time_elapsed = rstan::get_elapsed_time(model),
      divergent = rstan::get_num_divergent(model),
      bfmi = rstan::get_bfmi(model),
      ess = coda::effectiveSize(rstan::As.mcmc.list(model))[1:8],
      mv_rhat = coda::gelman.diag(rstan::As.mcmc.list(model))$mpsrf
    )
  })

  ## For old model
  oldConvergence <- lapply(oldFiles, function(i) {
    model <- readRDS(i)$fit
    list(
      time_elapsed = rstan::get_elapsed_time(model),
      divergent = rstan::get_num_divergent(model),
      bfmi = rstan::get_bfmi(model),
      ess = coda::effectiveSize(rstan::As.mcmc.list(model))[1:8],
      mv_rhat = coda::gelman.diag(rstan::As.mcmc.list(model))$mpsrf
    )
  })

  rawConvergence <- list(new = newConvergence,
                         old = oldConvergence)

  rm(newConvergence, oldConvergence)

  saveRDS(rawConvergence, file = "output/rawConvergence.RDS")
} else {
  rawConvergence <- readRDS(file = "output/rawConvergence.RDS")
}


# Time to converge --------------------------------------------------------

ttcData <- lapply(rawConvergence, function(m) {
  lapply(m, function(x) {
    as.data.table(x$time_elapsed, keep.rownames = "chain")
  }) |> rbindlist(idcol = "id")
}) |> rbindlist(idcol = "model")

ttcData[model == "new", file := gsub(pattern = "models/new/sim/|\\.RDS", replacement = "", newFiles[id])]
ttcData[model == "old", file := gsub(pattern = "models/old/sim/|\\.RDS", replacement = "", oldFiles[id])]

ttcData[, prior := fcase(
  grepl("narrow", file), "narrow",
  grepl("moderate", file), "moderate",
  grepl("wide", file), "wide"
)]

ttcData[, id := as.numeric(x = gsub("narrow|moderate|wide|Prior|\\-", "", file))]

ttcData[, prior := factor(prior,
                          levels = c("narrow","moderate","wide"),
                          labels = c("Narrow", "Moderate", "Wide"),
                          ordered = TRUE)]

ttcData[ , model := factor(model, levels = c("old", "new"), labels = c("Old", "New"))]

ttcData[, elapsed := warmup + sample]

p1 <- ggplot(ttcData, aes(elapsed * 1000, prior, fill = prior)) +
  facet_grid(rows = vars(model)) +
  tidybayes::stat_halfeye(adjust = 6, density = "unbounded", trim = FALSE, fill = "#069") +
  scale_y_discrete(expand = c(0.1,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(100,NA)) +
  labs(subtitle = "Time to Convergence (ms)", x = NULL, y = "Prior", fill = "Prior") +
  see::theme_modern(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.major.y = element_line())

ttcProp <- ttcData[, list(proportion = elapsed[model == "New"] / elapsed[model == "Old"]), list(prior, chain, id)]

p2 <- ggplot(ttcProp, aes(proportion, prior)) +
  tidybayes::stat_halfeye(adjust = 6, density = "unbounded", trim = F, fill = "#069") +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 1, col = "gray") +
  scale_y_discrete(expand = c(0.1,0)) +
  labs(subtitle = "Proportion of Time to Convergence (New/Old)", x = NULL, y = "Prior", fill = "Prior") +
  see::theme_modern(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.major.y = element_line())

fig_ttc <- ggpubr::ggarrange(p1, p2, nrow = 2, heights = c(2,1), common.legend = TRUE)

ggsave(filename = "output/ttc_fig.pdf", plot = fig_ttc, width = 6, height = 8)
