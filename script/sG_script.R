# Simulator sG

library(CurveRegistration)
library(protoABC)
library(parallel)
library(dplyr)

# Settings ----------------------

run_number = 1

set.seed(1)

Time <- seq(0, 200, by = 0.5)
alpha = seq(20, 180, by = 20)
theta = c(1, 0.7, 0.01)

alpha_norm   <- alpha[seq(1,length(alpha), by = 2)]
alpha_cauchy <- alpha[seq(2,length(alpha), by = 2)]

n_runs = 4000
pacc_final = 0.002

var_mat <- diag(c(9, 0.01^2))

y <- simulator_sGaussian(Time, param = theta, alpha_norm = alpha_norm, alpha_cauchy = alpha_cauchy)
y_kmmd <- EasyMMD::kmmd(y, var = var_mat)
x <- simulator_sGaussian(Time, param = theta, alpha_norm = alpha_norm, alpha_cauchy = alpha_cauchy)

distance_fun(y, x, registration = TRUE, distance = "MMD", y_kmmd = y_kmmd, var = var_mat, threshold = 6)

distance_args <- list(
  Time = Time,
  alpha_norm = alpha_norm,
  alpha_cauchy = alpha_cauchy,
  obs = y,
  var = var_mat,
  threshold = 6,
  registration = TRUE,
  distance = "MMD",
  method = "DP2",
  mean_global = 0,
  sigma_a = 2.5
)

distance_args$y_kmmd <- EasyMMD::kmmd(distance_args$obs, var = distance_args$var)

loss_sG(theta, distance_args)

# ABC -----------------

prior_sGaussian <- prior_unif(c(0, 0, 0), c(3, 2, 0.02), var_names = c("sigma[phi]", "rho[phi]", "sigma[epsilon]"), eval = FALSE)

prior_sGaussian_eval <- prior_unif(c(0, 0, 0), c(3, 2, 0.02), var_names = c("sigma[phi]", "rho[phi]", "sigma[epsilon]"), eval = TRUE)

cov_func <- function(x){
  robust::covRob(x)$cov
}

abc_control <- list(
  prior_eval = prior_sGaussian_eval,
  n = n_runs,
  pacc_final = pacc_final
)

#cl <- makeCluster(detectCores() - 1) #USER
cl <- "mclapply" #HPC

## sigma_a <- 5 -----------------

distance_args$sigma_a <- 2.5

### Registration <- FALSE -------------

distance_args$registration <- FALSE

#### Distance <- "FR" -------------------

distance_args$distance <- "FR"

ABC_sG_RF_FR <- abc_start(prior_sGaussian, loss_sG, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    sigma_a      = distance_args$sigma_a
  )

save(ABC_sG_RF_FR, file = "ABC_sG_RF_FR.RData")

summary(ABC_sG_RF_FR)

#### Distance <- "MMD" -------------------

distance_args$distance <- "MMD"

ABC_sG_RF_MD <- abc_start(prior_sGaussian, loss_sG, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    sigma_a      = distance_args$sigma_a
  )

save(ABC_sG_RF_MD, file = "ABC_sG_RF_MD.RData")

summary(ABC_sG_RF_MD)

### Registration <- TRUE -------------

distance_args$registration <- TRUE

#### Distance <- "FR" -------------------

distance_args$distance <- "FR"

ABC_sG_RT_FR <- abc_start(prior_sGaussian, loss_sG, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    sigma_a      = distance_args$sigma_a
  )

save(ABC_sG_RT_FR, file = "ABC_sG_RT_FR.RData")

summary(ABC_sG_RT_FR)

#### Distance <- "MMD" -------------------

distance_args$distance <- "MMD"

ABC_sG_RT_MD <- abc_start(prior_sGaussian, loss_sG, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    sigma_a      = distance_args$sigma_a
  )

save(ABC_sG_RT_MD, file = "ABC_sG_RT_MD.RData")

summary(ABC_sG_RT_MD)

# Packaging everything up

sessionInfo()

ABC_sG <- dplyr::bind_rows(ABC_sG_RF_FR, ABC_sG_RF_MD, ABC_sG_RT_FR, ABC_sG_RT_MD)

save(ABC_sG, file = "ABC_sG.RData")

save.image()




