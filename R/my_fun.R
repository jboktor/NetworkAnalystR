#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname my.fun
#' @export 
  my.fun <- function(){
    require('PCSF');
    edg <- get.edgelist(dat.in$data);
    edg <- as.data.frame(edg);
    edg$V3 <- rep(1, nrow(edg));
    colnames(edg) <- c("from", "to", "cost");
    ppi <- construct_interactome(edg);
    g <- PCSF(ppi, dat.in$terminals, w = 5, b = 100, mu = 0.0005);
    return(g);
  }
