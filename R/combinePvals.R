#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param pvalonesided PARAM_DESCRIPTION
#' @param nrep PARAM_DESCRIPTION
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
#' @rdname combinePvals
#' @export 
combinePvals <- function(pvalonesided,nrep,BHth=0.05, method) {
  listres=vector("list",3);
  nbstudies=length(pvalonesided);
  nbreptot=sum(nrep);
  if (nbreptot <2) {
    stop("Error: the argument \"nrep\" must be a vector with at least two values higher than 1")
  } 
  
  weight=sqrt(nrep/nbreptot);
  fstatmeta=function(g){
    vecptime=unlist(lapply(pvalonesided, FUN = function(x) x[g]));
    vec = qnorm(1 - vecptime);
    stattestg = sum(weight[1:length(pvalonesided)] * vec[1:length(pvalonesided)], na.rm = TRUE);
    stattestg;
  }
  
  fishersum <- function(pvec){
    return(sum(-2*log(pvec)))
  }
  
  if(method=="stouffer"){
    statpvalc=-unlist(lapply(rep(1:length(as.vector(pvalonesided[[1]])), 1), function(x) fstatmeta(x)));
    rpvalpvalc=2*(1-pnorm(abs(statpvalc)));
  }else{ # fisher
    data <- data.frame(pvalonesided);
    #data[data == 0] <- 1*10^-10;
    
    #note, p value are calculated for one side
    # pt (lower.tail=T by default) which tests if group A < group B
    # for one side
    fsum1 <- apply(data, 1, fishersum);
    rpvalpvalc1 = 1-pchisq(fsum1, df=(ncol(data)*2));
    
    # for the other side
    data <- 1-data;
    fsum2 <- apply(data, 1, fishersum);
    rpvalpvalc2 = 1-pchisq(fsum2, df=(ncol(data)*2));
    
    # report the min of two direction calculation
    rpvalpvalc <- pmin(rpvalpvalc1, rpvalpvalc2);
    
    # adding direction information
    statpvalc <- pmax(fsum1, fsum2);
    # if A<B sig, then it should be negative 
    statpvalc[statpvalc == fsum1]<- -statpvalc[statpvalc == fsum1];
  }
  
  rpvalpvalc <- p.adjust(rpvalpvalc,method="BH");
  res=which(rpvalpvalc<=BHth);
  listres[[1]]=res
  listres[[2]]=statpvalc;
  listres[[3]]=rpvalpvalc
  names(listres)=c("DEindices", "CombinedStat", "CombinedP")
  listres
}
