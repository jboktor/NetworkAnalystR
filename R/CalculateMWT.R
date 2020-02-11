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
#' @seealso 
#'  \code{\link[BiocGenerics]{unique}},\code{\link[BiocGenerics]{colSums}},\code{\link[BiocGenerics]{mean}},\code{\link[BiocGenerics]{cbind}}
#'  \code{\link[h2o]{Ops.H2OFrame}}
#'  \code{\link[stats]{TDist}}
#' @rdname CalculateMWT
#' @export 
#' @importFrom BiocGenerics unique rowMeans rowSums mean cbind
#' @importFrom h2o log
#' @importFrom stats pt
CalculateMWT <- function(xdat,grp,na.rm=TRUE){
  ## basic statistics
  glab = BiocGenerics::unique(grp)
  n1 = sum(grp==glab[1])
  n2 = sum(grp==glab[2])
  d1 = n1-1
  d2 = n2-1
  m1 = BiocGenerics::rowMeans(xdat[,grp==glab[1]], na.rm=na.rm)
  m2 = BiocGenerics::rowMeans(xdat[,grp==glab[2]], na.rm=na.rm)
  
  s2.g1 = BiocGenerics::rowSums((xdat[,grp==glab[1]]-m1)^2, na.rm=na.rm)/d1
  
  ## We might either have all NA in one group or variance = 0
  ## (e.g. might happen with RMA with small samples)
  ## In this situation we want to remove the gene
  s2.g1[s2.g1 == 0] <- NA
  
  s2.g2 = BiocGenerics::rowSums((xdat[,grp==glab[2]]-m2)^2, na.rm=na.rm)/d2
  s2.g2[s2.g2 == 0] <- NA
  
  ## If either s2.g1 or s2.g2 are NA this will be NA
  sig2 = (d1*s2.g1 + d2*s2.g2)/(d1+d2)
  fac = 1/n1 + 1/n2
  se2 = (sig2 * fac)
  
  ## F test
  
  lev.test = levene(xdat,grp)
  fFDR = lev.test$FDR
  fStat = lev.test$statistic
  
  
  ## ordinary Welch statistics
  se2.sep = s2.g1/n1 + s2.g2/n2
  df = se2.sep^2/((s2.g1/n1)^2/d1 + (s2.g2/n2)^2/d2)
  
  ## weighted formulas
  df.w  = fFDR*(d1+d2) + (1-fFDR)*df
  se2.w = fFDR*se2 + (1-fFDR)*se2.sep
  ds = est.hyper(z=h2o::log(se2.w),D=BiocGenerics::mean(df.w,na.rm=na.rm),d12=d1+d2)   
  
  ## ....................................... moderated Welch
  se2.com = (ds$d0*ds$s2 + df.w*se2.w)/(ds$d0 + df.w)
  Wm.t = (m1-m2)/sqrt(se2.com) ## Welch t
  df.com = ds$d0 + df.w      ## df
  
  
  Wm.pval = stats::pt(-abs(Wm.t), df= df.com) * 2
  
  ## ................. Compute Global FDR
  
  fdr <- NULL
  
  return(list(MWT= Wm.t, coefficients=BiocGenerics::cbind((m1-m2)),pvalue = Wm.pval))
  
}
