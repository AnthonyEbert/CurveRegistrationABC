

#' @export
distance_fun <- function(obs, sim, registration, distance, method = "DP2", y_kmmd = NULL, var = NULL, threshold = Inf){

  # registration == 1: elastic distance
  # registration == 2: MMD
  # registration == 3: R-MMD

  if(registration){
    reg_output <- register_functions(obs, sim, method = method)
    sim[,2] <- reg_output$value
  }

  if(distance == "FR"){
    output <- sqrt(pracma::trapz(obs[,1], (fdasrvf::f_to_srvf(obs[,2], obs[,1]) - fdasrvf::f_to_srvf(sim[,2], sim[,1]))^2))
  }

  if(distance == "MMD"){
    output <- EasyMMD::MMD(obs, sim, y_kmmd = y_kmmd, bias = TRUE, var = var, threshold = threshold)
  }

  return(output)
}
