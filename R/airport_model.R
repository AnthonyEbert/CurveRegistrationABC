

#' @export
sim_airport <- function(params, inp){
  global_input <- AirportSim::global_fun_1(mu = params[1], vm2 = params[2])

  nat_input <- AirportSim::airport_list_1$nat_level
  nat_input$rate_imm <- params[3:4]

  flight_input <- inp$flight_level

  if(inp$flight_effect){
    flight_input$arrive <- runif(length(flight_input$arrive), pmax.int(0, flight_input$arrive -20), flight_input$arrive + 20)
  }

  input <- AirportSim::airport_fun_1(global_level = global_input, nat_level = nat_input, flight_level = flight_input)

  Passenger_df <- do.call(AirportSim::AirportSimulate1, input)

  xy_obs_syn <- AirportSim::post_process_1(Passenger_df)

  x_full <- AirportSim::convert_stamps_to_hist(xy_obs_syn, breaks = inp$breaks)

  names(x_full) <- c("key", "Time", "y")

  x_ac <- x_full[which(x_full[,1] == "x_obs"),]
  x_imm <- x_full[which(x_full[,1] == "y_obs"),]

  x_ac <- as.matrix(x_ac[,-1])
  x_imm <- as.matrix(x_imm[,-1])

  output <- list(x_ac = x_ac, x_imm = x_imm, Passenger_df = Passenger_df)

  return(output)
}

#' @export
loss_airport <- function(params, inp){

  obs <- inp$obs
  obs_ac <- obs$x_ac
  obs_imm <- obs$x_imm

  sim <- sim_airport(params, inp)
  sim_ac <- sim$x_ac
  sim_imm <- sim$x_imm

  if(inp$registration | inp$correction){
    reg_output_ac <- register_functions(obs_ac, sim_ac, method = inp$method, warp_ret = TRUE)
  }

  if(inp$correction){

    sim_imm <- correction_airport(sim$Passenger_df, reg_output_ac$warp_fun, inp$breaks)

  }

  if(inp$registration){
    sim_ac[,2] <- reg_output_ac$value
  }

  output_ac <- distance_fun(obs_ac, sim_ac, registration = FALSE, distance = inp$distance, method = inp$method, y_kmmd = NULL, var = inp$var, threshold = inp$threshold)

  output_imm <- distance_fun(obs_imm, sim_imm, registration = inp$registration_imm, distance = inp$distance, method = inp$method, y_kmmd = NULL, var = inp$var, threshold = inp$threshold)

  return(as.numeric(output_ac + output_imm))
}


#' @export
correction_airport <- function(sim_Passenger_df, warp_fun, breaks){

  sim_out <- sim_Passenger_df %>%
    group_by(nat) %>%
    mutate(
      depart_imm = queuecomputer::queue(warp_fun(arrive_imm), service_imm, server_imm[[1]])
    )

  sim_out <- sim_out %>%
    AirportSim::post_process_1() %>%
    AirportSim::convert_stamps_to_hist(breaks = breaks)

  sim_imm <- sim_out %>%
    filter(key == "y_obs") %>%
    select(-key) %>%
    as.matrix()

  return(sim_imm)
}









