# note, entrez.vec could contain NA/null, cannot use rownames
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
#' @rdname doEntrezIDAnot
#' @export 
doEntrezIDAnot <- function(entrez.vec){
  gene.map <-  queryGeneDB("entrez", data.org);
  
  hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
  anot.mat <- gene.map[hit.inx, c("gene_id", "symbol", "name")];
  
  na.inx <- is.na(hit.inx);
  anot.mat[na.inx, "symbol"] <- entrez.vec[na.inx];
  anot.mat[na.inx, "name"] <- 'NA';
  return(anot.mat);
}
