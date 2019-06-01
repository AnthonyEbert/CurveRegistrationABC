
library(airGR)
library(protoABC)
library(parallel)
library(CurveRegistration)
library(dplyr)

# Settings -------------------


data(L0123001)

InputsModel <-
  CreateInputsModel(
    FUN_MOD = RunModel_GR4J,
    DatesR = BasinObs$DatesR,
    Precip = BasinObs$P,
    PotEvap = BasinObs$E
  )

Ind_Run <-
  seq(which(format(BasinObs$DatesR, format = "%d/%m/%Y") == "01/03/1997"),
      which(format(BasinObs$DatesR, format = "%d/%m/%Y") == "01/01/1998"))

RunOptions <- CreateRunOptions(
  FUN_MOD = RunModel_GR4J,
  InputsModel = InputsModel,
  IndPeriod_Run = Ind_Run,
  IniStates = NULL,
  IniResLevels = NULL,
  IndPeriod_WarmUp = NULL
)

Param <- c(257, 1, 88, 10, 0.05)

#Param <- c(115, 1.07, 128, 2.8)

OutputsModel <-
  RunModel_GR4J(InputsModel = InputsModel,
                RunOptions = RunOptions,
                Param = Param[1:4])

var_mat <- diag(c(100, 1e-2))

obs_syn <- simulator_hydro(Param, InputsModel, RunOptions)

distance_args <- list(
  InputsModel = InputsModel,
  RunOptions = RunOptions,
  obs = obs_syn,
  var = var_mat,
  y_kmmd = EasyMMD::kmmd(obs_syn, var = var_mat),
  registration = FALSE,
  distance = "FR",
  threshold = 6,
  method = "DP2"
)

loss_hydro(Param, distance_args)

n_runs <- 10000
pacc_final <- 0.02

# ABC -----------------

prior_hydro <- protoABC::prior_unif(lhs_lim = c(100, -5, 20, 1.1, 0.01), rhs_lim = c(1200, 3, 300, 30, 0.08), var_names = c("x1", "x2", "x3", "x4", "sigma"))

prior_eval_hydro <- protoABC::prior_unif(lhs_lim = c(100, -5, 20, 1.1, 0.01), rhs_lim = c(1200, 3, 300, 30, 0.08), var_names = c("x1", "x2", "x3", "x4", "sigma"), eval = TRUE)

abc_control <- list(
  n = n_runs,
  prior_eval = prior_eval_hydro,
  pacc_final = pacc_final
)

#cl <- makeCluster(detectCores() - 1) #USER
cl <- "mclapply" #HPC

## Synthetic data ------------------

### Registration <- FALSE -------------

distance_args$registration <- FALSE

#### Distance <- "FR" ------------------

distance_args$distance <- "FR"

ABC_hydro_syn_RF_FR <- protoABC::abc_start(prior_hydro, loss_hydro, distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    data = "synthetic",
    registration = distance_args$registration,
    distance = distance_args$distance
  )

save(ABC_hydro_syn_RF_FR, file = "ABC_hydro_syn_RF_FR.RData")

summary(ABC_hydro_syn_RF_FR)

#### Distance <- "MMD" ------------------

distance_args$distance <- "MMD"

ABC_hydro_syn_RF_MD <- protoABC::abc_start(prior_hydro, loss_hydro, distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    data = "synthetic",
    registration = distance_args$registration,
    distance = distance_args$distance
  )

save(ABC_hydro_syn_RF_MD, file = "ABC_hydro_syn_RF_MD.RData")

summary(ABC_hydro_syn_RF_MD)

### Registration <- TRUE -----------

distance_args$registration <- TRUE

#### Distance <- "FR" ------------------

distance_args$distance <- "FR"

ABC_hydro_syn_RT_FR <- protoABC::abc_start(prior_hydro, loss_hydro, distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    data = "synthetic",
    registration = distance_args$registration,
    distance = distance_args$distance
  )

save(ABC_hydro_syn_RT_FR, file = "ABC_hydro_syn_RT_FR.RData")

summary(ABC_hydro_syn_RT_FR)

#### Distance <- "MMD" ------------------

distance_args$distance <- "MMD"

ABC_hydro_syn_RT_MD <- protoABC::abc_start(prior_hydro, loss_hydro, distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    data = "synthetic",
    registration = distance_args$registration,
    distance = distance_args$distance
  )

save(ABC_hydro_syn_RT_MD, file = "ABC_hydro_syn_RT_MD.RData")

summary(ABC_hydro_syn_RT_MD)

## Observed data -----------

n_days <- length(Ind_Run)
day_num <- seq(1, n_days)

obs_true <- as.matrix(data.frame(day_num = day_num, value = BasinObs$Qmm[Ind_Run]))

distance_args$obs <- obs_true
distance_args$y_kmmd <- EasyMMD::kmmd(obs_true, var = var_mat)

loss_hydro(Param, distance_args)

### Registration <- FALSE -------------

distance_args$registration <- FALSE

#### Distance <- "FR" ------------------

distance_args$distance <- "FR"

ABC_hydro_obs_RF_FR <- protoABC::abc_start(prior_hydro, loss_hydro, distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    data = "synthetic",
    registration = distance_args$registration,
    distance = distance_args$distance
  )

save(ABC_hydro_obs_RF_FR, file = "ABC_hydro_obs_RF_FR.RData")

summary(ABC_hydro_obs_RF_FR)

#### Distance <- "MMD" ------------------

distance_args$distance <- "MMD"

ABC_hydro_obs_RF_MD <- protoABC::abc_start(prior_hydro, loss_hydro, distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    data = "synthetic",
    registration = distance_args$registration,
    distance = distance_args$distance
  )

save(ABC_hydro_obs_RF_MD, file = "ABC_hydro_obs_RF_MD.RData")

summary(ABC_hydro_obs_RF_MD)

### Registration <- TRUE -----------

distance_args$registration <- TRUE

#### Distance <- "FR" ------------------

distance_args$distance <- "FR"

ABC_hydro_obs_RT_FR <- protoABC::abc_start(prior_hydro, loss_hydro, distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    data = "synthetic",
    registration = distance_args$registration,
    distance = distance_args$distance
  )

save(ABC_hydro_obs_RT_FR, file = "ABC_hydro_obs_RT_FR.RData")

summary(ABC_hydro_obs_RT_FR)

#### Distance <- "MMD" ------------------

distance_args$distance <- "MMD"

ABC_hydro_obs_RT_MD <- protoABC::abc_start(prior_hydro, loss_hydro, distance_args, method = "RABC", control = abc_control, cl = cl) %>%
  mutate(
    data = "synthetic",
    registration = distance_args$registration,
    distance = distance_args$distance
  )

save(ABC_hydro_obs_RT_MD, file = "ABC_hydro_obs_RT_MD.RData")

summary(ABC_hydro_obs_RT_MD)


sessionInfo()

ABC_hydro_syn <- dplyr::bind_rows(ABC_hydro_syn_RF_FR, ABC_hydro_syn_RF_MD, ABC_hydro_syn_RT_FR, ABC_hydro_syn_RT_MD)

ABC_hydro_obs <- dplyr::bind_rows(ABC_hydro_obs_RF_FR, ABC_hydro_obs_RF_MD, ABC_hydro_obs_RT_FR, ABC_hydro_obs_RT_MD)

save(ABC_hydro_syn, file = "ABC_hydro_syn.RData")
save(ABC_hydro_obs, file = "ABC_hydro_obs.RData")

save.image(file = "ABC_hydro_everything.RData")


