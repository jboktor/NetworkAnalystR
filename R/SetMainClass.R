# note, setup the main class, keep the original order
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cls.lbl PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetMainClass
#' @export 
SetMainClass<-function(cls.lbl){
  lbls <- as.character(dataSet$meta.info[[cls.lbl]]);
  lvls.orig <- unique(lbls);
  cls <- factor(lbls, levels=lvls.orig, ordered=T);
  dataSet$cls <- cls; # record main cls
  dataSet <<- dataSet;
  return(levels(cls));
}
