
library(AirportSim)
library(dplyr)
library(CurveRegistration)
library(protoABC)

cov_func <- function(x){
  robust::covRob(x)$cov
}

set.seed(3)

flight_level <- AirportSim::generate_flightlevel(5, 1000)

true_params <- c(0.02, 0.64, 0.25, 0.4)

distance_args <- list(
  flight_effect = TRUE,
  registration = FALSE,
  registration_imm = FALSE,
  var = diag(c(2, 2)),
  method = "DP",
  threshold = 6,
  breaks = seq(0, 1000, by = 5),
  flight_level = flight_level,
  correction = FALSE,
  distance = "FR"
)

output <- sim_airport(true_params, distance_args)

distance_args$obs <- output

loss_airport(true_params, distance_args)

# ABC -----------

param_names <- c("mu", "vm2", "lambda_f", "lambda_l")

prior_airport <- protoABC::prior_unif(c(0, 0, 0, 0), c(0.05, 1, 1, 1), var_names = param_names, eval = FALSE)

prior_eval_airport <- protoABC::prior_unif(c(0, 0, 0, 0), c(0.05, 1, 1, 1), var_names = param_names, eval = TRUE)

cov_func <- function(x){
  robust::covRob(x)$cov
}

abc_control <- list(
  n          = 4000,
  pacc_final = 0.02,
  prior_eval = prior_eval_airport,
  cov_func   = cov_func
)

#cl <- parallel::makeCluster(parallel::detectCores() - 1)
#parallel::clusterEvalQ(cl = cl, expr = library(dplyr))

cl <- "mclapply"

## Registration <- FALSE

distance_args$registration <- FALSE
distance_args$registration_imm <- FALSE

#### Distance <- "FR" -------------------

distance_args$distance <- "FR"

##### correction <- FALSE

distance_args$correction <- FALSE

ABC_airport_RF_FR_CF <- abc_start(prior_airport, loss_airport, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    correction   = distance_args$correction
  )

save(ABC_airport_RF_FR_CF, file = "ABC_airport_RF_FR_CF.RData")

##### correction <- TRUE

distance_args$correction <- TRUE

ABC_airport_RF_FR_CT <- abc_start(prior_airport, loss_airport, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    correction   = distance_args$correction
  )

save(ABC_airport_RF_FR_CT, file = "ABC_airport_RF_FR_CT.RData")

#### Distance <- "FR" -------------------

distance_args$distance <- "MMD"

##### correction <- FALSE

distance_args$correction <- FALSE

ABC_airport_RF_MD_CF <- abc_start(prior_airport, loss_airport, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    correction   = distance_args$correction
  )

save(ABC_airport_RF_MD_CF, file = "ABC_airport_RF_MD_CF.RData")

##### correction <- TRUE

distance_args$correction <- TRUE

ABC_airport_RF_MD_CT <- abc_start(prior_airport, loss_airport, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    correction   = distance_args$correction
  )

save(ABC_airport_RF_MD_CT, file = "ABC_airport_RF_MD_CT.RData")

## Registration <- TRUE

distance_args$registration <- TRUE
distance_args$registration_imm <- FALSE

#### Distance <- "FR" -------------------

distance_args$distance <- "FR"

##### correction <- FALSE

distance_args$correction <- FALSE

ABC_airport_RT_FR_CF <- abc_start(prior_airport, loss_airport, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    correction   = distance_args$correction
  )

save(ABC_airport_RT_FR_CF, file = "ABC_airport_RT_FR_CF.RData")

##### correction <- TRUE

distance_args$correction <- TRUE

ABC_airport_RT_FR_CT <- abc_start(prior_airport, loss_airport, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    correction   = distance_args$correction
  )

save(ABC_airport_RT_FR_CT, file = "ABC_airport_RT_FR_CT.RData")

#### Distance <- "MMD" -------------------

distance_args$distance <- "MMD"

##### correction <- FALSE

distance_args$correction <- FALSE

ABC_airport_RT_MD_CF <- abc_start(prior_airport, loss_airport, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    correction   = distance_args$correction
  )

save(ABC_airport_RT_MD_CF, file = "ABC_airport_RT_MD_CF.RData")

##### correction <- TRUE

distance_args$correction <- TRUE

ABC_airport_RT_MD_CT <- abc_start(prior_airport, loss_airport, distance_args = distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    registration = distance_args$registration,
    distance     = distance_args$distance,
    correction   = distance_args$correction
  )

save(ABC_airport_RT_MD_CT, file = "ABC_airport_RT_MD_CT.RData")

sessionInfo()

ABC_airport <- dplyr::bind_rows(ABC_airport_RF_FR_CF, ABC_airport_RF_FR_CT, ABC_airport_RF_MD_CF, ABC_airport_RF_MD_CT, ABC_airport_RT_FR_CF, ABC_airport_RT_FR_CT, ABC_airport_RT_MD_CF, ABC_airport_RT_MD_CT)

save(ABC_airport, file = "ABC_airport.RData")
save.image("ABC_airport_everything.RData")

