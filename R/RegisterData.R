# When multiple genelists/datasets/results, record their name and save the data as .RDS file
# a) Current dataSet object
# Note, the memory will only contain one dataSet object. By default, the last one will be the current dataSet object;
# Users can switch this (from the interface) to specify which data is load into memory (dataSet object)
# b) Include for certain analysis
# For chord and heatmap analysis, users can do multiple selection (include)
# All datasets are selected by default (1 for selected, 0 for unselected)
# note, dataSet need to have "name" property
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataSet PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname RegisterData
#' @export 
RegisterData <- function(dataSet){
  dataName <- dataSet$name;
  saveRDS(dataSet, file=dataName);
  if(!is.null(dataSet$data.raw)){ # save memory for meta-analysis mode
    dataSet$data.raw <- NULL;
  }
  dataSet <<- dataSet;
  if(anal.type == "metadata"){
    mdata.all[[dataName]] <<- 1;
  }else{
    mdata.all <<- lapply(mdata.all, function(x){ x <- 0;});
    mdata.all[[dataName]] <<- 1;
  }
  return(1);
}
