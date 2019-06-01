
#' @export
simulator_hydro <- function(param, InputsModel, RunOptions, one_d = FALSE, lambda_bc = 0.5, A_bc = 5){

  OutputsModel <-
    airGR::RunModel_GR4J(InputsModel = InputsModel,
                         RunOptions = RunOptions,
                         Param = param[1:4])

  n_days <- length(OutputsModel$Qsim)

  day_num <- seq(1, n_days)

  tvalue = box_cox(OutputsModel$Qsim, lambda_bc, A_bc) + rnorm(n_days, 0, param[5])

  output <- as.matrix(data.frame(day_num = day_num, value = inv_box_cox(tvalue, lambda_bc, A_bc)))

  if(one_d){
    output <- output[,2]
  }

  return(output)
}

#' @export
loss_hydro <- function(param, inp){

  sim <- simulator_hydro(param, inp$InputsModel, inp$RunOptions)

  output <- distance_fun(inp$obs, sim, registration = inp$reg, distance = inp$distance, method = inp$method,  y_kmmd = inp$y_kmmd, var = inp$var, threshold = inp$threshold)

  return(output)
}




#' @export
box_cox <- function(x, lambda, A){

  output <- ((x + A)^lambda - 1)/lambda

  return(output)
}

#' @export
inv_box_cox <- function(x, lambda, A){

  output <- (x * lambda + 1)^(1/lambda) - A

  return(output)
}






