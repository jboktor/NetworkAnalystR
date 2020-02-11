#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param id.vec PARAM_DESCRIPTION
#' @param idType PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doAnnotation
#' @export 
doAnnotation <- function(id.vec, idType){
  if(idType %in% c("entrez", "symbol", "refseq", "genbank", "emblgene","emblprotein", "embltranscript", "orfid", "tair", "wormbase")){
    anot.id <- doGeneIDMapping(id.vec, idType);
  }else{
    anot.id <- doProbeMapping(id.vec, idType);
    names(anot.id) <- id.vec;
  }
  return(anot.id);        
}
