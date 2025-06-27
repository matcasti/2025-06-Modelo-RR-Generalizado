library(data.table)

#> ## Data generation process

#> - Describe how synthetic RRi data were generated for the simulations. Specify that data were created using the original model formulation with a predefined set of "true" underlying physiological parameters.
#> - Detail the range and distribution of these "true" parameters, ensuring they represent a variety of plausible physiological scenarios.

## Function to generate normal samples from a range, where the limits are on the 99% CI limit
normal_from_range <- function(n, min, max, theoretical = FALSE) {
  if (theoretical) {
    c(mean = (max-min)/2 + min, sd = (max-min)/(2 * qnorm(.99)))
  } else {
    rnorm(n, mean = (max-min)/2 + min, sd = (max-min)/(2 * qnorm(.99)))
  }
}

original_model <- function(t, alpha, beta, c, lambda, phi, tau, delta) {
  alpha -
    beta / (1 + exp(-lambda * (t - tau))) +
    (c * beta) / (1 + exp(-phi * (t - tau - delta)))
}

do <- function(.dt, ...) {
  expr <- substitute(.dt[...])
  eval(expr)
}

n <- 100

## Draw true parameters from Uniform distribution
set.seed(1234)
params <- data.table(
  id = seq_len(n),
  alpha = normal_from_range(n, 600, 1200),
  beta = normal_from_range(n, 100, 500),
  c = normal_from_range(n, 0.5, 1.1),
  lambda = normal_from_range(n, 1, 7),
  phi = normal_from_range(n, 1, 7),
  tau = normal_from_range(n, 3, 8),
  delta = normal_from_range(n, 3, 9),
  sigma = normal_from_range(n, 20, 150)
)

list(
  alpha = normal_from_range(n, 600, 1200, theoretical = TRUE),
  beta = normal_from_range(n, 100, 500, theoretical = TRUE),
  c = normal_from_range(n, 0.5, 1.1, theoretical = TRUE),
  lambda = normal_from_range(n, 1, 7, theoretical = TRUE),
  phi = normal_from_range(n, 1, 7, theoretical = TRUE),
  tau = normal_from_range(n, 3, 8, theoretical = TRUE),
  delta = normal_from_range(n, 3, 9, theoretical = TRUE),
  sigma = normal_from_range(n, 20, 150, theoretical = TRUE)
)

sim_data <- params |>
  do(, {
    t <- seq(0.1, 15, 0.1)
    rri_true <- original_model(t, alpha, beta, c, lambda, phi, tau, delta)
    rri_obs <- round(rri_true + rnorm(length(rri_true), 0, sigma))
    bar_rri <- mean(rri_obs)
    s_rri <- sd(rri_obs)
    list(rri_true, rri_obs, t, bar_rri, s_rri)
  }, id)

library(ggplot2)

ggplot(sim_data, aes(x = t, colour = ordered(id))) +
  geom_line(aes(y = rri_obs), show.legend = FALSE) +
  #scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(title = "Observed RRi Curve", subtitle = "Based on Simulated Parameters",
       y = "RRi (ms)", x = "Time (min)") +
  theme_classic()

#> - Crucially, explain how **varying time scales** were incorporated into the simulated data. Provide specific examples of different total protocol durations (e.g., 5 min, 15 min, 60 min, 120 min) and/or varying sampling rates (e.g., 1 Hz, 10 Hz) to simulate diverse experimental designs.

library(brms)

model_formula <- bf(rri_obs ~ alpha -
                      beta / (1 + exp(-lambda * (t - tau))) +
                      (c * beta) / (1 + exp(-phi * (t - tau - delta))),
                    alpha + beta + c + lambda + phi+  tau + delta ~ 1,
                    nl = TRUE)

model_prior <- c(
  set_prior("normal(900, 128.9575)", nlpar = "alpha", class = "b"),
  set_prior("normal(300, 85.97166)", nlpar = "beta", class = "b"),
  set_prior("normal(0.8, 0.1289575)", nlpar = "c", class = "b"),
  set_prior("normal(4, 1.289575)", nlpar = "lambda", class = "b"),
  set_prior("normal(4, 1.289575)", nlpar = "phi", class = "b"),
  set_prior("normal(5.5, 1.074646)", nlpar = "tau", class = "b"),
  set_prior("normal(5.5, 1.289575)", nlpar = "delta", class = "b"),
  set_prior("normal(85, 27.94079)", class = "sigma", lb = 0)
)


models <- vector("list", length = n)
times <- numeric(length = n)
for (i in seq_len(n)) {
  model_data <- sim_data[id == i]

  times[i] <- system.time(
    models[[i]] <- brm(
      formula = model_formula,
      data = model_data,
      family = gaussian(),
      prior = model_prior,
      iter = 10000, warmup = 5000,
      chains = 4, cores = 4
    )
  )
}



#> - Describe the method of **noise introduction** to mimic real-world data (e.g., addition of Gaussian noise with a specified standard deviation, or noise characteristics derived from empirical RRi data).



#> - State the number of Monte Carlo simulations performed for each unique combination of parameters, time scale, and noise level (e.g., $N=1000$ simulations per scenario).
