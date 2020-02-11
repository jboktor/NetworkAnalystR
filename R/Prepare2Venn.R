#3
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dat PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Prepare2Venn
#' @export 
Prepare2Venn <- function(dat){
  nms <- names(dat);  
  a <- nms[1];
  b <- nms[2];
  ab <- paste(a, b, sep="");
  
  a.l <- dat[[a]];
  b.l <- dat[[b]];
  
  vennData <- list();
  vennData[[a]] <- setdiff(a.l, b.l);
  vennData[[b]] <- setdiff(b.l, a.l);    
  vennData[[ab]] <- intersect(b.l, a.l);
  vennData <<- vennData;
}
