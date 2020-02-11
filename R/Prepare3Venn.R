#7
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
#' @rdname Prepare3Venn
#' @export 
Prepare3Venn <- function(dat){
  nms <- names(dat);  
  a <- nms[1];
  b <- nms[2];
  c <- nms[3];
  ab <- paste(a, b, sep="");
  ac <- paste(a, c, sep="");
  bc <- paste(b, c, sep="");
  abc <- paste(a, b, c, sep="");
  
  a.l <- dat[[a]];
  b.l <- dat[[b]];
  c.l <- dat[[c]];
  
  vennData <- list();
  vennData[[a]] <- setdiff(a.l, union(b.l, c.l));
  vennData[[b]] <- setdiff(b.l, union(a.l, c.l));    
  vennData[[c]] <- setdiff(c.l, union(a.l, b.l));    
  vennData[[ab]] <- setdiff(intersect(a.l, b.l), c.l);
  vennData[[ac]] <- setdiff(intersect(a.l, c.l), b.l);
  vennData[[bc]] <- setdiff(intersect(b.l, c.l), a.l);
  vennData[[abc]] <- intersect(intersect(a.l, b.l), c.l);
  vennData <<- vennData;
}
