#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param xdat PARAM_DESCRIPTION
#' @param grp PARAM_DESCRIPTION
#' @param na.rm PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname levene
#' @export 
levene <- function (xdat, grp, na.rm = TRUE) 
{
  glab = unique(grp)
  ngr = length(glab)
  n = mn = s2 = NULL
  X0 = NULL
  for (i in 1:ngr) {
    ndx = grp == glab[i]
    mni = rowMeans(xdat[, ndx], na.rm = na.rm)
    x0 = xdat[, ndx] - mni
    X0 = cbind(X0, x0)
  }
  xdat = abs(X0)
  for (i in 1:ngr) {
    ndx = grp == glab[i]
    ni = sum(ndx)
    mni = rowMeans(xdat[, ndx], na.rm = na.rm)
    x0 = xdat[, ndx] - mni
    s2i = rowSums(x0 * x0, na.rm = na.rm)/(ni - 1)
    n = c(n, ni)
    mn = cbind(mn, mni)
    s2 = cbind(s2, s2i)
  }
  N = sum(n)
  mmn = rowSums(xdat, na.rm = na.rm)/N
  mn0 = mn - mmn
  num = rowSums(t(t(mn0 * mn0) * n))/(ngr - 1)
  den = rowSums(t(t(s2 * (n - 1))))/(N - ngr)
  F3 = num/den
  pval = pf(F3, df1 = ngr - 1, df2 = N - ngr)
  pval = ifelse(pval < 0.5, 2 * pval, 2 * (1 - pval))
  lvn.FDR = pval2FDR(pval)
  return(list(statistic = F3, pvalue = pval, FDR = lvn.FDR))
}
