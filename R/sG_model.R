
#' @export
loss_sG <- function(param, inp){

  obs <- inp$obs

  if(!isTRUE(inp$skew)){
    sim <- simulator_sGaussian(inp$Time, param = param, alpha_norm = inp$alpha_norm, alpha_cauchy = inp$alpha_cauchy, mean_global = inp$mean_global, sigma_a = inp$sigma_a)
  } else {
    sim <- simulator_skew(inp$Time, param = param, alpha_norm = inp$alpha_norm, mean_global = inp$mean_global, sigma_a = inp$sigma_a)
  }

  out <- distance_fun(obs, sim, registration = inp$registration, distance = inp$distance, method = inp$method, y_kmmd = inp$y_kmmd, var = inp$var, threshold = 6)

  return(as.numeric(out))
}

#' Simulated example 1 (Gaussian and Cauchy peaks)
#' @examples
#' Time <- seq(0, 200, by = 0.5)
#' alpha = seq(20, 180, by = 20)
#' theta = c(1, 0.7, 0.01)
#' alpha_norm   <- alpha[seq(1,length(alpha), by = 2)]
#' alpha_cauchy <- alpha[seq(2,length(alpha), by = 2)]
#' y <- simulator_sGaussian(Time, param = theta, alpha_norm = alpha_norm, alpha_cauchy = alpha_cauchy)
#' @export
simulator_sGaussian <- function(Time, param = c(1, 0.7, 0.01), alpha_norm = c(20, 60), alpha_cauchy = c(40, 80), mean_global = 0, sigma_a = 5){

  sigma_global  <- param[1]
  scale_global  <- param[2]
  sigma_e       <- param[3]

  mean_global   <- ifelse(is.na(param[4]), mean_global, param[4])
  sigma_a       <- ifelse(is.na(param[5]), sigma_a, param[5])

  a_norm        <- rnorm(length(alpha_norm), mean_global, sigma_a)
  a_cauchy      <- rnorm(length(alpha_cauchy), mean_global, sigma_a)

  y <- rep(NA, length(Time))

  for(i in 1:length(Time)){
    y[i] <- sum(dnorm(Time[i], mean = a_norm + alpha_norm, sd = sigma_global)) +
      sum(dcauchy(Time[i], location = a_cauchy + alpha_cauchy, scale = scale_global)) +
      rnorm(1, 0, sigma_e)
  }

  output <- as.matrix(data.frame(Time = Time, Output = y))

  return(output)
}

#' Simulated example 2 (Gaussian skewed peaks)
#' @examples
#' Time <- seq(0, 200, by = 0.5)
#' alpha = seq(20, 180, by = 20)
#' theta = c(1, 0.7, 1)
#' y <- simulator_skew(Time, alpha_norm = alpha, param = theta)
#' @export
simulator_skew <- function(Time, param = c(1, 0.7, 0.01), alpha_norm = c(20, 40, 60, 80), mean_global = 0, sigma_a = 5){

  sigma_global  <- param[1]
  skew_global   <- param[2]
  sigma_e       <- param[3]

  mean_global   <- ifelse(is.na(param[4]), mean_global, param[4])
  sigma_a       <- ifelse(is.na(param[5]), sigma_a, param[5])

  a_norm        <- rnorm(length(alpha_norm), mean_global, sigma_a)

  y <- matrix(nrow = length(Time), ncol = length(alpha_norm))

  for(i in 1:length(Time)){
    for(j in 1:length(alpha_norm)){
      y[i,j] <- sn::dsn(Time[i], dp = sn::cp2dp(cp = c(alpha_norm[j] + a_norm[j], param[1:2]), "SN")) + rnorm(1, 0, sigma_e)
    }
  }

  y <- rowSums(y)

  output <- as.matrix(data.frame(Time = Time, Output = y))

  return(output)
}









