#combine effect size
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ES PARAM_DESCRIPTION
#' @param varES PARAM_DESCRIPTION
#' @param BHth PARAM_DESCRIPTION, Default: 0.05
#' @param method PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname combineES
#' @export 
combineES <- function (ES, varES, BHth = 0.05, method){
  if(method == "rem"){
    useREM = TRUE;
  }else{
    useREM = FALSE;
  }
  
  num.studies <- dim(ES)[2];
  
  Qvals <- f.Q.NA(ES, varES)
  if (useREM) {
    varES <- varES + tau2.NA(Qvals, num.studies, my.weights = 1/varES)
  }
  wt <- 1/varES
  MUvals <- rowSums(ES * wt, na.rm = TRUE)/rowSums(wt, na.rm = TRUE)
  MUsES <- sqrt(abs(1/rowSums(wt, na.rm = TRUE)))
  zSco <- MUvals/MUsES
  rpvalESc = 2 * (1 - pnorm(abs(zSco)))
  res = which(p.adjust(rpvalESc, method = "BH") <= BHth);
  listres <- list();
  listres[[1]] = res
  listres[[2]] = zSco
  listres[[3]] = MUvals; # pool effect size
  listres[[4]] = wt; # wt for each studies, it is matrix with one column for each studies
  names(listres) = c("DEindices", "TestStatistic", "PooledEffectSize", "Weights")
  listres
}
