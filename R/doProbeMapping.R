# from probe ID to entrez ID 
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param probe.vec PARAM_DESCRIPTION
#' @param platform PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doProbeMapping
#' @export 
doProbeMapping <- function(probe.vec, platform){
  platform.path <- paste(lib.path,  data.org, "/", platform, ".rds", sep="");
  probe.map <- readRDS(platform.path);
  if(is.null(probe.vec)){
    entrez <- probe.map[, "entrez"];
  }else{
    hit.inx <- match(probe.vec, probe.map[, "probe"]);
    entrez <- probe.map[hit.inx, "entrez"];
  }
  return(entrez);
}
