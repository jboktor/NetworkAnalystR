#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param tstat PARAM_DESCRIPTION
#' @param ntilde PARAM_DESCRIPTION
#' @param m PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname effectsize
#' @export 
effectsize <- function(tstat,ntilde,m){
  cm=gamma(m/2)/(sqrt(m/2)*gamma((m-1)/2))
  d=tstat/sqrt(ntilde)
  dprime=cm*d
  terme1=m/((m-2)*ntilde)
  vard=terme1+d^2*(terme1*ntilde-1/cm^2)
  vardprime=cm^2*(terme1+dprime^2*(terme1*ntilde-1/cm^2))
  result=cbind(d,vard,dprime,vardprime)
  colnames(result)=c("d","vard","dprime","vardprime")
  result
}
