

# Model implementations ---------------------------------------------------

orig_model <- function(t, params, ...) {
  with(params, {
    alpha -
      beta / (1 + exp(-lambda * (t - tau))) +
      (c * beta) / (1 + exp(-phi * (t - tau - delta)))
  })
}

new_model <- function(t, params, ...) {
  with(params, {
    fun_alpha(alphaR, rriDelta, rriMin) -
      fun_beta(betaR, rriDelta) / (1 + exp(-fun_lambda(lambdaR) * (t - fun_tau(tauR, tDelta, tMin)))) +
      (fun_c(cR) * fun_beta(betaR, rriDelta)) / (1 + exp(-fun_phi(phiR) * (t - fun_tau(tauR, tDelta, tMin) - fun_delta(deltaR, tDelta))))
  })
}

# Parameter functions -----------------------------------------------------

## From transformed to original parameters
fun_alpha <- function(alphaR, rriDelta, rriMin) {
  rriMin + (2 * rriDelta) * (exp(alphaR) / (1 + exp(alphaR)))
}
fun_beta <- function(betaR, rriDelta) {
  (2 * rriDelta) * (exp(betaR) / (1 + exp(betaR)))
}
fun_c <- function(cR) {
  (2 * exp(cR)) / (1 + exp(cR))
}
fun_lambda <- function(lambdaR) {
  exp(lambdaR)
}
fun_phi <- function(phiR) {
  exp(phiR)
}
fun_tau <- function(tauR, tDelta, tMin) {
  tMin + tDelta * (exp(tauR) / (1 + exp(tauR)))
}
fun_delta <- function(deltaR, tDelta) {
  tDelta * (exp(deltaR) / (1 + exp(deltaR)))
}

## From original to transformed parameters
fun_alphaR <- function(alpha, rriDelta, rriMin) {
  log((alpha - rriMin)/(2 * rriDelta - alpha + rriMin))
}
fun_betaR <- function(beta,  rriDelta) {
  log(beta/(2 * rriDelta - beta))
}
fun_cR <- function(c) {
  log(c / (2 - c))
}
fun_lambdaR <- function(lambda) {
  log(lambda)
}
fun_phiR <- function(phi) {
  log(phi)
}
fun_tauR <- function(tau, tDelta, tMin) {
  log((tau - tMin)/(tDelta - tau + tMin))
}
fun_deltaR <- function(delta, tDelta) {
  log(delta/(tDelta - delta))
}

# Functions to estimate transformed `sigma` -------------------------------

sigma_alphaR <- function(sigma_alpha, alpha, rriDelta, rriMin) {
  deriv <- (1/(2 * rriDelta - alpha + rriMin) + (alpha - rriMin) / (2 * rriDelta - alpha + rriMin)^2) /
    ((alpha - rriMin) / (2 * rriDelta - alpha + rriMin))
  abs(deriv) * sigma_alpha
}

sigma_betaR <- function(sigma_beta, beta, rriDelta) {
  deriv <- (1/(2 * rriDelta - beta) + beta / (2 * rriDelta - beta)^2) /
    (beta/(2 * rriDelta - beta))
  abs(deriv) * sigma_beta
}

sigma_cR <- function(sigma_c, c) {
  abs((1/(2 - c) + c/(2 - c)^2)/(c/(2 - c))) * sigma_c
}

sigma_lambdaR <- function(sigma_lambda, lambda) {
  abs(1 / lambda) * sigma_lambda
}

sigma_phiR <- function(sigma_phi, phi) {
  abs(1 / phi) * sigma_phi
}

sigma_tauR <- function(sigma_tau, tau, tDelta, tMin) {
  deriv <- (1 / (tDelta - tau + tMin) + (tau - tMin) / (tDelta - tau + tMin)^2) /
    ((tau - tMin) / (tDelta - tau + tMin))
  abs(deriv) * sigma_tau
}

sigma_deltaR <- function(sigma_delta, delta, tDelta) {
  deriv <- (1 / (tDelta - delta) + delta / (tDelta - delta)^2) /
    (delta/(tDelta - delta))
  abs(deriv) * sigma_delta
}


# Auxiliary functions ----------------------------------------------------

do <- function(.dt, ...) {
  expr <- substitute(.dt[...])
  eval(expr)
}

# Data generation functions -----------------------------------------------

normal <- function(n, distParams = list(location = 0, scale = 1)) {
  rnorm(n = n,
        mean = distParams$location,
        sd = distParams$scale)
}

uniform <- function(n, distParams = list(location = 0, scale = 1), mult = 2) {
  runif(n = n,
        min = distParams$location - distParams$scale * mult,
        max = distParams$location + distParams$scale * mult)
}


# Prior generation functions ----------------------------------------------

modelPrior <- function(parameter, mult = 1) {
  paste0("normal(",parameter$location,", ",parameter$scale * mult,")")
}
