##################################################
## R scripts for NetworkAnalyst 
## Description: Meta Analysis Methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################
# for multiple class, only select two
# also record all grp lbls
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataName PARAM_DESCRIPTION
#' @param grps PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetGroupContrast
#' @export 
SetGroupContrast <- function(dataName, grps){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  if(length(levels(dataSet$cls))>2){ 
    print("Updating group contrasts .....");
    grp.nms <- strsplit(grps, " vs. ")[[1]];
    sel.inx <- as.character(dataSet$cls) %in% grp.nms;
    
    # regenerate factor to drop levels, force the levels order
    group <- factor(dataSet$cls[sel.inx], levels=grp.nms);  
    data <- dataSet$data[, sel.inx];
    dataSet$cls <- group;
    dataSet$data <- data;
    RegisterData(dataSet);  
  }
}
