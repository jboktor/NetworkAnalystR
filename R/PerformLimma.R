# perfor differential analysis for array/RNA seq data
# for two groups only (used for meta-analysis)
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param group PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformLimma
#' @export 
PerformLimma<-function(data, group){
  require(limma);
  data <- data;
  design <- model.matrix(~-1 + group);
  fit = lmFit(data, design)
  
  grps.cmp <- paste("group", levels(group)[2], " - ", "group", levels(group)[1], sep="");
  myargs <- list(grps.cmp, levels = design);
  contrast.matrix <- do.call(makeContrasts, myargs);
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit);
  gc();
  return (list(fit.obj=fit));
}
