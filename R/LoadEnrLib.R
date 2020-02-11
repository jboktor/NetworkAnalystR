##################################################
## R script for NetworkAnalyst
## Description: Gene/Probe/Protein ID Annotation
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fun.type PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname LoadEnrLib
#' @export 
LoadEnrLib <- function(fun.type){
  if(tolower(fun.type) == 'kegg'){ 
    current.geneset <- LoadKEGGLib();
  }else if(tolower(fun.type) == 'reactome'){ 
    current.geneset <-LoadREACTOMELib();
  }else if(tolower(fun.type) == 'motif'){ 
    current.geneset <- LoadMotifLib();
  }else{ # GO
    current.geneset <-LoadGOLib(fun.type);
  }
  return(current.geneset)
}
