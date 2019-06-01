#' @export
register_functions <- function(obs, sim, method, warp_ret = FALSE){

  reg_output <- fdasrvf::pair_align_functions(obs[,2], sim[,2], time = obs[,1], method = method)

  values <- reg_output$f2tilde

  warp_fun <- NULL

  if(warp_ret){
    warp_fun <- warper(sim[,1], reg_output$gam)
  }

  return(list(values = values, gam = reg_output$gam, warp_fun = warp_fun))

}
