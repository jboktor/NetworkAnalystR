#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param deMethod PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetupDesignMatrix
#' @export 
SetupDesignMatrix<-function(deMethod){
  cls <- dataSet$cls; 
  design <- model.matrix(~ 0 + cls) # no intercept
  colnames(design) <- levels(cls);
  dataSet$design <- design;
  dataSet$de.method <- deMethod;
  dataSet <<- dataSet;
  
  return(1);
}
