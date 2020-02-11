##################################################
## R script for NetworkAnalyst
## Description: functions only for list data analysis
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fileName PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SaveListHeatmapJSON
#' @export 
SaveListHeatmapJSON <- function(fileName){
  if(numOfLists>1){
    json.res <- PrepareMultiListHeatmapJSON();
  }else{
    json.res <- PrepareListHeatmapJSON();
  }
  require(RJSONIO);
  json.mat <- toJSON(json.res, .na='null');
  sink(fileName);
  cat(json.mat);
  sink();
  current.msg <<- "Data is now ready for heatmap visualization!";
  return(1);
}
