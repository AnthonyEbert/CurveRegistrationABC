#' @export
warper <- function(x, gam){
  ord = order(x)
  t = x[ord]
  gam = gam[ord]

  return(
    function(t_eval){

      gam2 <- t[1] + gam * (t[length(t)] - t[1])
      x <- approx(gam2, t, t_eval)$y
      x[which(is.na(x))] <- t_eval[which(is.na(x))]
      return(x)
    })
}
