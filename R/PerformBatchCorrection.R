# this is a transient R environment for specific function, then close
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param method PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformBatchCorrection
#' @export 
PerformBatchCorrection <- function(method){
  
  print("performing combat batch correction ....");
  
  library(RSclient);
  rsc <- RS.connect();
  
  RS.assign(rsc, "my.dir", getwd()); 
  RS.eval(rsc, setwd(my.dir));
  
  my.fun <- function(){
    library('sva');
    inmex.meta <- readRDS("inmex_meta.rds");
    data.lbl <- inmex.meta$data.lbl
    pheno <- data.frame(cbind(inmex.meta$cls.lbl, data.lbl));
    modcombat <- model.matrix(~1, data=pheno)
    batch <- data.lbl;
    inmex.meta$data <- ComBat(dat=inmex.meta$data, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE);
    saveRDS(inmex.meta, "inmex_meta.rds");
    return(1);
  }
  RS.assign(rsc, my.fun);
  res <-  RS.eval(rsc, my.fun());
  RS.close(rsc);
  
  return(res);
}
