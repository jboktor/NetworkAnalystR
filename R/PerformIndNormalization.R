# note, here also update data type array/count
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataName PARAM_DESCRIPTION
#' @param norm.opt PARAM_DESCRIPTION
#' @param auto.opt PARAM_DESCRIPTION
#' @param dataType PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformIndNormalization
#' @export 
PerformIndNormalization <- function(dataName, norm.opt, auto.opt, dataType){
  
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  msg <- NULL;
  data <- dataSet$data.orig;
  data <- PerformDataNormalization(data, norm.opt);
  if(length(data)==1 && data == 0){
    return(0);
  }
  msg <- paste(norm.msg, msg);
  
  if(auto.opt==1){
    row.nms <- rownames(data);
    col.nms <- colnames(data);
    data<-apply(data, 2, AutoNorm);
    msg <- paste(msg, "Autoscaling performed.", collapse=" ");
    rownames(data) <- row.nms;
    colnames(data) <- col.nms;
  }
  
  dataSet$data <- data;
  dataSet$type <- dataType;
  RegisterData(dataSet);
  current.msg <<- msg;
  return(1);
}
