# Simulator skew

library(CurveRegistration)
library(protoABC)
library(parallel)
library(dplyr)

# Settings ----------------------

run_number = 1

set.seed(1)

Time <- seq(0, 200, by = 0.5)
alpha = seq(20, 180, by = 20)
theta = c(1, 0.75, 0.01)

n_runs = 1000
pacc_final = 0.01

var_mat <- diag(c(9, 0.01^2))

y <- simulator_skew(Time, param = theta, alpha_norm = alpha, sigma_a = 2.5)
y_kmmd <- EasyMMD::kmmd(y, var = var_mat)
x <- simulator_skew(Time, param = theta, alpha_norm = alpha, sigma_a = 2.5)

distance_fun(y, x, registration = TRUE, distance = "MMD", y_kmmd = y_kmmd, var = var_mat, threshold = 6)

distance_args <- list(
  Time = Time,
  alpha_norm = alpha,
  obs = y,
  var = var_mat,
  threshold = 6,
  registration = TRUE,
  distance = "MMD",
  method = "DP2",
  mean_global = 0,
  sigma_a = 2.5,
  skew = TRUE
)

distance_args$y_kmmd <- EasyMMD::kmmd(distance_args$obs, var = distance_args$var)

loss_sG(theta, distance_args)

# ABC -----------------

prior_skew <- prior_unif(c(0, -0.9, 0), c(3, 0.9, 0.02), var_names = c("sigma[phi]", "eta", "sigma[epsilon]"), eval = FALSE)

prior_skew_eval <- prior_unif(c(0, -0.9, 0), c(3, 0.9, 0.02), var_names = c("sigma[phi]", "eta", "sigma[epsilon]"), eval = TRUE)

cov_func <- function(x){
  robust::covRob(x)$cov
}

abc_control <- list(
  prior_eval = prior_skew_eval,
  n = n_runs,
  pacc_final = pacc_final,
  a = 0.5
)

#cl <- makeCluster(detectCores() - 1) #USER
cl <- "mclapply" #HPC

### Registration <- FALSE -------------

distance_args$registration <- FALSE

#### Distance <- "FR" -------------------

distance_args$distance <- "FR"

ABC_skew_RF_FR <- abc_start(prior_skew, loss_sG, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    sigma_a      = distance_args$sigma_a
  )

save(ABC_skew_RF_FR, file = "ABC_skew_RF_FR.RData")

summary(ABC_skew_RF_FR)

#### Distance <- "MMD" -------------------

distance_args$distance <- "MMD"

ABC_skew_RF_MD <- abc_start(prior_skew, loss_sG, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    sigma_a      = distance_args$sigma_a
  )

save(ABC_skew_RF_MD, file = "ABC_skew_RF_MD.RData")

summary(ABC_skew_RF_MD)

### Registration <- TRUE -------------

distance_args$registration <- TRUE

#### Distance <- "FR" -------------------

distance_args$distance <- "FR"

ABC_skew_RT_FR <- abc_start(prior_skew, loss_sG, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    sigma_a      = distance_args$sigma_a
  )

save(ABC_skew_RT_FR, file = "ABC_skew_RT_FR.RData")

summary(ABC_skew_RT_FR)

#### Distance <- "MMD" -------------------

distance_args$distance <- "MMD"

ABC_skew_RT_MD <- abc_start(prior_skew, loss_sG, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    sigma_a      = distance_args$sigma_a
  )

save(ABC_skew_RT_MD, file = "ABC_skew_RT_MD.RData")

summary(ABC_skew_RT_MD)

# Packaging everything up

sessionInfo()

ABC_skew <- dplyr::bind_rows(ABC_skew_RF_FR, ABC_skew_RF_MD, ABC_skew_RT_FR, ABC_skew_RT_MD)

save(ABC_skew, file = "ABC_skew.RData")

save.image()




