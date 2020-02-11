#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param pval PARAM_DESCRIPTION
#' @param lim PARAM_DESCRIPTION, Default: 0.7
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname pval2FDR
#' @export 
pval2FDR <-function (pval, lim = 0.7) 
{
  n1 = length(pval)
  ok.id <- 1:n1
  if (any(is.na(pval))) {
    ok.id <- which(!is.na(pval))
    pval <- na.omit(pval)
  }
  n = length(pval)
  Fp = rank(pval)/length(pval)
  p0 = sum(pval > lim)/((1 - lim) * n)
  p0 = min(p0, 1)
  FDRp = p0 * pmin(pval/Fp, 1)
  ord = order(pval)
  FDR.o = FDRp[ord]
  b = rev(cummin(rev(FDR.o)))
  FDR = rep(0, n)
  FDR[ord] = b
  out.FDR <- rep(NA, n1)
  out.FDR[ok.id] <- FDR
  attr(out.FDR, "p0") <- p0
  return(out.FDR)
}
