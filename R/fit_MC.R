
# Manila Clam Source Code

#' Condition an operating model for Manila Clam
#'
#' Uses the RCM model of OpenMSE to fit an operating model to data
#'
#' @param In A list object of class ('In') that includes the operating model parameters (slot 2, class OM) and data (slot 3, class RCMinput) for conditioning
#' @param sims Integer or vector of integers - the number of simulations for the operating model (e.g. 96) or the specific vector of simulations for the operating model (e.g. 13:48)
#' @param max_F Positive real number - the maximum instantaneous mortality rate in any historical time step
#' @param comp_like The likelihood function used for composition data (age / length) c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2") see ?RCM
#' @param resample Logical - should parameters be drawn from the variance covariance matrix of the MLE fit (T) or just the MLE parameter values (F)
#' @param parallel Logical - should fitting use parallel processing
#' @param silent Logical - should all RCM messages be repressed (ie not read out onto the console)?
#' @examples
#' cond.MC("In.MC.E")
#' @author T. Carruthers
#' @seealso \link{RCM}
#' @export
cond.MC = function(RCMinput, sims = 12, max_F = 2.3, comp_like = "multinomial", resample = F, parallel = F,silent=T){
  cores = 1
  if(parallel){
    setup()
    cores=parallel::detectCores()/2
  }
  if(length(sims) == 1) simy = 1:sims
  if(length(sims) > 1) simy = sims
  OM = SubCpars(RCMinput$OM, simy)
  RCMfit = RCM(OM, RCMinput[[3]], selectivity = "logistic_length", s_selectivity=c("logistic_age","B","B"),
             max_F = max_F, mean_fit = T, comp_like = comp_like, condition = "catch", cores = cores,
             drop_nonconv=T,drop_highF=T,resample=resample,LWT=list(Index=c(1/8,1,1)), silent=silent)

  RCMfit@OM@Name = paste("Manila_Clam CMA", RCMinput[[1]])
  RCMfit

}
