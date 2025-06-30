
# Prepare workspace -------------------------------------------------------

## Load libraries
library(data.table)
library(brms)
library(ggplot2)

## Load custom functions
original_model <- function(t, alpha, beta, c, lambda, phi, tau, delta) {
  alpha -
    beta / (1 + exp(-lambda * (t - tau))) +
    (c * beta) / (1 + exp(-phi * (t - tau - delta)))
}

do <- function(.dt, ...) {
  expr <- substitute(.dt[...])
  eval(expr)
}

normal <- function(n, distParams = list(location = 0, scale = 1)) {
  rnorm(n, distParams$location, distParams$scale)
}

uniform <- function(n, distParams = list(location = 0, scale = 1)) {
  runif(n, distParams$location - distParams$scale, distParams$location + distParams$scale)
}


# Create synthetic data ---------------------------------------------------

n <- 10 ## Number of Monte Carlo simulations

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
  delta = list(location = 3.5, scale = 1),
  sigma = list(location = 30, scale = 5)
)

set.seed(1234) ## Seed for reproducibility

## Draw true underlying parameters from normal distribution
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
    sigma = normal(n, sigma)
  )
})

## Generate synthetic true and observed (with gaussian noise) RRi signal
sim_data <- simParams |>
  do(, {
    t <- seq(0.1, 15, 0.01)
    rri_true <- original_model(t, alpha, beta, c, lambda, phi, tau, delta)
    rri_obs <- round(rri_true + rnorm(n = length(rri_true), sd = sigma))
    bar_rri <- mean(rri_obs)
    s_rri <- sd(rri_obs)
    list(rri_true, rri_obs, t, bar_rri, s_rri)
  }, id)

## Visualize simulated RRi data
ggplot(sim_data, aes(x = t, colour = ordered(id))) +
  geom_line(aes(y = rri_obs), show.legend = FALSE) +
  scale_x_continuous(expand = c(0,0)) +
  labs(title = "Observed RRi Signal", subtitle = "From Simulated Parameters",
       y = "RRi (ms)", x = "Time (min)") +
  theme_classic()



model_formula <- bf(rri_obs ~ alpha -
                      beta / (1 + exp(-lambda * (t - tau))) +
                      (c * beta) / (1 + exp(-phi * (t - tau - delta))),
                    alpha + beta + c + lambda + phi+  tau + delta ~ 1,
                    nl = TRUE)

priorParams <- lapply(trueParams, function(x) {
  with(x, paste0("normal(",location,", ",scale,")"))
})

model_prior <- c(
  set_prior(priorParams$alpha, nlpar = "alpha", class = "b", lb = 0),
  set_prior(priorParams$beta, nlpar = "beta", class = "b", lb = 0),
  set_prior(priorParams$c, nlpar = "c", class = "b", lb = 0),
  set_prior(priorParams$lambda, nlpar = "lambda", class = "b", lb = 0),
  set_prior(priorParams$phi, nlpar = "phi", class = "b", lb = 0),
  set_prior(priorParams$tau, nlpar = "tau", class = "b", lb = 0),
  set_prior(priorParams$delta, nlpar = "delta", class = "b", lb = 0),
  set_prior(priorParams$sigma, class = "sigma", lb = 0)
)

times <- vector("list", length = n)
for (i in seq_len(n)) {
  model_data <- sim_data[id == i]

  gc()
  times[[i]] <- system.time(
    brm(
      formula = model_formula,
      data = model_data,
      family = gaussian(),
      prior = model_prior,
      iter = 10000, warmup = 5000,
      chains = 4, cores = 4,
      save_model = "models/sim-original/model.stan",
      file = paste0("models/sim-original/sim-",i,".RDS"),
      a
    )
  )
}


