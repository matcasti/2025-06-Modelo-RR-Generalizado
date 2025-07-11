
# Script info -------------------------------------------------------------

## File: 1-sim.R
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

## Load custom functions
source("R/_functions.R")

## Decide if start without any
remove_models <- FALSE
remove_precompiled <- FALSE

if (remove_models) {
  status <- file.remove(
    list.files(path = c("models/new/sim","models/old/sim"),
               pattern = "\\.RDS$",
               full.names = TRUE,
               recursive = TRUE)
  )
  if (any(!status)) {
    warning("Some models were not deleted!")
  } else if (length(status) == 0) {
    message("All models were already deleted!")
  }
  rm(status)
}

if (remove_precompiled) {
  status <- file.remove(
    list.files(path = c("models/new","models/old"),
               pattern = "Prior\\.RDS$",
               full.names = TRUE,
               recursive = TRUE)
  )
  if (any(!status)) {
    warning("Some prior-only models were not deleted!")
  } else if (length(status) == 0) {
    message("All prior-only models were already deleted!")
  }
  rm(status)
}

# 1.1. Estimate parameters from the literature ------------------------------

## Establish true parameters based on Castillo-Aguilar et al. (2025):
## - Location based on rounded estimates from validation paper.
## - Scale based on 10 times the SE from validated estimates.
trueParams <- list(
  alpha = list(location = 850, scale = 50),
  beta = list(location = 350, scale = 70),
  c = list(location = 0.85, scale = 0.1),
  lambda = list(location = 3, scale = 0.6),
  phi = list(location = 2.5, scale = 0.6),
  tau = list(location = 6.5, scale = 0.5),
  delta = list(location = 3.5, scale = 0.5),
  sigma = list(location = 30, scale = 5)
)

# 1.2. Generate Synthetic RRi -----------------------------------------------

n <- 1000 ## Number of Monte Carlo simulations
## Note: the actual number of models fitted will be n x 6, given that for each
## n there are two models to test (old and new), and three priors to test.

t <- seq(0.01, 15, 0.01) ## Time vector

set.seed(1234) ## Seed for reproducibility

## Draw true underlying true parameters from normal distribution
## using true distributional parameters for each model parameter
simParams <- with(trueParams, {
  data.table(
    id = seq_len(n),
    alpha = normal(n, alpha),
    beta = normal(n, beta),
    c = normal(n, c),
    lambda = normal(n, lambda),
    phi = normal(n, phi),
    tau = normal(n, tau),
    delta = normal(n, delta),
    sigma = normal(n, sigma),
    tMin = min(t),
    tDelta = diff(x = range(t))
  )
})

## Generate synthetic true and observed (with Gaussian noise) RRi signal
simData <- simParams |>
  do(j = {
    rri_true <- orig_model(t, list(alpha = alpha, beta = beta, c = c,
                                  lambda = lambda, phi = phi, tau = tau,
                                  delta = delta))

    rri_obs <- rri_true + rnorm(n = length(rri_true), sd = sigma)

    ## Estimate rriDelta and rriMin from observed data
    rriDelta <- diff(x = range(rri_obs))
    rriMin <- min(rri_obs)

    list(rri_true, rri_obs, t, rriDelta, rriMin, tMin, tDelta)
  }, by = id)

# 2.1. Compute new model parameters from old model parameters ---------------

## Match subject-wise the rriDelta and rriMin from the simulation data (simData).
## Note: These parameters are needed to compute alphaR and betaR
simParams <- simParams[simData[, .SD[1], id][, list(id, rriDelta, rriMin)], on = "id"]

simParams[, `:=`(
  alphaR = fun_alphaR(alpha, rriDelta, rriMin),
  betaR = fun_betaR(beta, rriDelta),
  cR = fun_cR(c),
  lambdaR = fun_lambdaR(lambda),
  phiR = fun_phiR(phi),
  tauR = fun_tauR(tau, tDelta, tMin),
  deltaR = fun_deltaR(delta, tDelta)
), id]

# 2.2. Estimate true location and scale for each new parameter --------------

trueParams <- within(trueParams, {
  rriDelta <- mean(simParams$rriDelta)
  rriMin <- mean(simParams$rriMin)

  tMin <- unique(simParams$tMin)
  tDelta <- unique(simParams$tDelta)

  alphaR <- list(location = fun_alphaR(alpha$location, rriDelta, rriMin),
                 scale = sigma_alphaR(alpha$scale, alpha$location, rriDelta, rriMin))
  betaR <- list(location = fun_betaR(beta$location, rriDelta),
                scale = sigma_betaR(beta$scale, beta$location, rriDelta))
  cR <- list(location = fun_cR(c$location),
             scale = sigma_cR(c$scale, c = c$location))
  lambdaR <- list(location = fun_lambdaR(lambda$location),
                  scale = sigma_lambdaR(lambda$scale, lambda = lambda$location))
  phiR <- list(location = fun_phiR(phi$location),
               scale = sigma_phiR(phi$scale, phi = phi$location))
  tauR <- list(location = fun_tauR(tau$location, tDelta, tMin),
               scale = sigma_tauR(tau$scale, tau = tau$location, tDelta, tMin))
  deltaR <- list(location = fun_deltaR(delta$location, tDelta),
                 scale = sigma_deltaR(delta$scale, delta = delta$location, tDelta))
})

saveRDS(trueParams, file = "output/trueParams.RDS")

# -------------------------------------------------------------------------

if (FALSE) {
  ## Visualize simulated RRi data, only if R is running interactively
  ggplot2::ggplot(simData, ggplot2::aes(x = t, colour = ordered(id))) +
    ggplot2::geom_line(ggplot2::aes(y = rri_obs), show.legend = FALSE) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::labs(title = "Observed RRi Signal", subtitle = "From Simulated Parameters",
         y = "RRi (ms)", x = "Time (min)") +
    ggplot2::theme_classic()
}

# -------------------------------------------------------------------------

# 3.1. Generate narrow, moderate and wide priors ----------------------------
# Rationale: to test convergence and model fitness at varying prior-information
# knowledge

## Multiplier for model standard deviation as model priors
priorMult <- list(narrow = 1/2, moderate = 1, wide = 2)

## For each multiplier, define model parameters' priors:

## For the new model parameters
newPrior <- lapply(priorMult, function(x) {
  with(trueParams, {
    set_prior(modelPrior(alphaR, mult = x), nlpar = "alphaR", class = "b") +
      set_prior(modelPrior(betaR, mult = x), nlpar = "betaR", class = "b") +
      set_prior(modelPrior(cR, mult = x), nlpar = "cR", class = "b") +
      set_prior(modelPrior(lambdaR, mult = x), nlpar = "lambdaR", class = "b") +
      set_prior(modelPrior(phiR, mult = x), nlpar = "phiR", class = "b") +
      set_prior(modelPrior(tauR, mult = x), nlpar = "tauR", class = "b") +
      set_prior(modelPrior(deltaR, mult = x), nlpar = "deltaR", class = "b") +
      set_prior(modelPrior(sigma, mult = x), class = "sigma", lb = 0)
  })
})

saveRDS(newPrior, file = "models/new/newPrior.RDS")

## For the old model parameters
oldPrior <- lapply(priorMult, function(x) {
  with(trueParams, {
    set_prior(modelPrior(alpha, mult = x), nlpar = "alpha", class = "b", lb = 0) +
      set_prior(modelPrior(beta, mult = x), nlpar = "beta", class = "b", lb = 0) +
      set_prior(modelPrior(c, mult = x), nlpar = "c", class = "b", lb = 0) +
      set_prior(modelPrior(lambda, mult = x), nlpar = "lambda", class = "b", lb = 0) +
      set_prior(modelPrior(phi, mult = x), nlpar = "phi", class = "b", lb = 0) +
      set_prior(modelPrior(tau, mult = x), nlpar = "tau", class = "b", lb = 0) +
      set_prior(modelPrior(delta, mult = x), nlpar = "delta", class = "b", lb = 0) +
      set_prior(modelPrior(sigma, mult = x), class = "sigma", lb = 0)
  })
})

saveRDS(oldPrior, file = "models/old/oldPrior.RDS")

# 4.1. Define new and old model formulas ------------------------------------

newFormula <- bf(
  formula = rri_obs ~ (rriMin + 2 * rriDelta * (exp(alphaR) / (1 + exp(alphaR)))) -
    (2 * rriDelta * (exp(betaR) / (1 + exp(betaR)))) / (1 + exp(-exp(lambdaR) * (t - (tMin + tDelta * (exp(tauR) / (1 + exp(tauR))))))) +
    (((2 * exp(cR)) / (1 + exp(cR))) * (2 * rriDelta * (exp(betaR) / (1 + exp(betaR))))) / (1 + exp(-exp(phiR) * (t - (tMin + tDelta * (exp(tauR) / (1 + exp(tauR)))) - (tDelta * (exp(deltaR) / (1 + exp(deltaR))))))),

  alphaR + betaR + cR + lambdaR + phiR + tauR + deltaR ~ 1,
  nl = TRUE
)

saveRDS(newFormula, file = "models/new/newFormula.RDS")

oldFormula <- bf(
  formula = rri_obs ~ alpha -
    beta / (1 + exp(-lambda * (t - tau))) +
    (c * beta) / (1 + exp(-phi * (t - tau - delta))),

  alpha + beta + c + lambda + phi + tau + delta ~ 1,
  nl = TRUE
)

saveRDS(oldFormula, file = "models/old/oldFormula.RDS")

# 4.2. Fit prior-only models  -----------------------------------------------
# Rationale: to avoid to recompile the model at each iteration

## For each model prior...
priorTypes <- c(narrow = "narrow", moderate = "moderate", wide = "wide")

model_newPrior <- lapply(priorTypes, function(x) {
  brm(
    formula = newFormula,
    data = simData[id == 1],
    family = gaussian(),
    prior = newPrior[[x]],
    iter = 10000, warmup = 5000,
    chains = 4, cores = 4, seed = 1234,
    sample_prior = "only",
    file = paste0("models/new/only_",x,"Prior.RDS")
  )
})

## Fit a prior-only model to avoid recompiling the model at each time
model_oldPrior <- lapply(priorTypes, function(x) {
  brm(
    formula = oldFormula,
    data = simData[id == 1],
    family = gaussian(),
    prior = oldPrior[[x]],
    iter = 10000, warmup = 5000,
    chains = 4, cores = 4, seed = 1234,
    sample_prior = "only",
    file = paste0("models/old/only_",x,"Prior.RDS")
  )
})


# -------------------------------------------------------------------------

if (FALSE) {
  conditional_effects(model_oldPrior$moderate, effects = "t")
  conditional_effects(model_newPrior$moderate, effects = "t")
}

# -------------------------------------------------------------------------

# 5.1. Fit models to data for each subject and for each prior ---------------

## For each model prior...
priorTypes

## And for each individual
subjectId <- unique(simData$id)

gc(full = TRUE)

message("Now starting simulation on model: NEW")

## New model run
model_newPosterior <- lapply(priorTypes, function(i) {
  gc(full = TRUE)
  message("Now starting with prior: ", i)
  lapply(subjectId, function(j) {
    message("Now starting with subject: ", j)
    mod_file <- paste0("models/new/sim/",formatC(j, digits = 4, width = 4, flag = "0"),"-",i,"Prior.RDS")
    if (!file.exists(mod_file)) {
      tryCatch(
        expr = {
          brm(
            formula = newFormula,
            data = simData[id == j],
            family = gaussian(),
            prior = newPrior[[i]],
            fit = model_newPrior[[i]],
            iter = 10000, warmup = 5000,
            chains = 4, cores = 4, seed = 1234,
            file = mod_file
          ); "Success"
        },
        error = function(e) {
          as.character(e)
        },
        warning = function(w) {
          as.character(w)
        }
      )
    } else {
      "Ommited"
    }
  })
})

unlist(model_newPosterior) |> table()

gc(full = TRUE)

message("Now starting simulation on model: OLD")

## Old model run
model_oldPosterior <- lapply(priorTypes, function(i) {
  gc(full = TRUE)
  message("Now starting with prior: ", i)
  lapply(subjectId, function(j) {
    message("Now starting with subject: ", j)
    mod_file <- paste0("models/old/sim/",formatC(j, digits = 4, width = 4, flag = "0"),"-",i,"Prior.RDS")
    if (!file.exists(mod_file)) {
      tryCatch(
        expr = {
          brm(
            formula = oldFormula,
            data = simData[id == j],
            family = gaussian(),
            prior = oldPrior[[i]],
            fit = model_oldPrior[[i]],
            iter = 10000, warmup = 5000,
            chains = 4, cores = 4, seed = 1234,
            file = mod_file
          ); "Success"
        },
        error = function(e) {
          as.character(e)
        },
        warning = function(w) {
          as.character(w)
        }
      )
    } else {
      "Ommited"
    }
  })
})

unlist(model_oldPosterior) |> table()

message("All models fitted!")


# End of script -----------------------------------------------------------
