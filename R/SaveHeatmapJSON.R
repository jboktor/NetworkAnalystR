# retrun the json obj
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
#' @rdname SaveHeatmapJSON
#' @export 
SaveHeatmapJSON <- function(fileName){
  
  if(anal.type == "metadata"){
    json.res <- PrepareMetaHeatmapJSON();
  }else{
    json.res <- PrepareExpressHeatmapJSON();
  }
  
  require(RJSONIO);
  json.mat <- toJSON(json.res, .na='null');
  sink(fileName);
  cat(json.mat);
  sink();
  current.msg <<- "Data is now ready for heatmap visualization!";
  return(1);
}
