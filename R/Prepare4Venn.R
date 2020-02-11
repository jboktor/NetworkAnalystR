# 15
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
#' @rdname Prepare4Venn
#' @export 
Prepare4Venn <- function(dat){
  nms <- names(dat);  
  a <- nms[1];
  b <- nms[2];
  c <- nms[3];
  d <- nms[4];
  ab <- paste(a, b, sep="");
  ac <- paste(a, c, sep="");
  ad <- paste(a, d, sep="");
  bc <- paste(b, c, sep="");
  bd <- paste(b, d, sep="");
  cd <- paste(c, d, sep="");
  abc <- paste(a, b, c, sep="");
  abd <- paste(a, b, d, sep="");
  acd <- paste(a, c, d, sep="");
  bcd <- paste(b, c, d, sep="");
  abcd <- paste(a, b, c, d, sep="");
  
  a.l <- dat[[a]];
  b.l <- dat[[b]];
  c.l <- dat[[c]];
  d.l <- dat[[d]];
  
  vennData <- list();
  vennData[[a]] <- setdiff(a.l, unique(c(b.l, c.l, d.l)));
  vennData[[b]] <- setdiff(b.l, unique(c(a.l, c.l, d.l)));    
  vennData[[c]] <- setdiff(c.l, unique(c(a.l, b.l, d.l)));    
  vennData[[d]] <- setdiff(d.l, unique(c(a.l, b.l, c.l))); 
  vennData[[ab]] <- setdiff(intersect(a.l, b.l), union(c.l, d.l));
  vennData[[ac]] <- setdiff(intersect(a.l, c.l), union(b.l, d.l));
  vennData[[ad]] <- setdiff(intersect(a.l, d.l), union(b.l, c.l));
  vennData[[bc]] <- setdiff(intersect(b.l, c.l), union(a.l, d.l));
  vennData[[bd]] <- setdiff(intersect(b.l, d.l), union(a.l, c.l));
  vennData[[cd]] <- setdiff(intersect(c.l, d.l), union(a.l, b.l));
  vennData[[abc]] <- setdiff(intersect(intersect(a.l, b.l), c.l), d.l);
  vennData[[abd]] <- setdiff(intersect(intersect(a.l, b.l), d.l), c.l);
  vennData[[acd]] <- setdiff(intersect(intersect(a.l, c.l), d.l), b.l);
  vennData[[bcd]] <- setdiff(intersect(intersect(b.l, c.l), d.l), a.l);
  vennData[[abcd]] <- intersect(intersect(a.l, b.l), intersect(c.l, d.l));
  vennData <<- vennData;
}
