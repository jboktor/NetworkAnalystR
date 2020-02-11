#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param entrez.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doID2LabelMapping
#' @export 
doID2LabelMapping <- function(entrez.vec){
  if(exists("nodeListu")){
    hit.inx <- match(entrez.vec, nodeListu[, "Id"]);
    symbols <- nodeListu[hit.inx, "Label"];
    
    # if not gene symbol, use id by itself
    na.inx <- is.na(symbols);
    symbols[na.inx] <- entrez.vec[na.inx];
    return(symbols);
  }else{ # network upload
    return(entrez.vec);
  }
} 
