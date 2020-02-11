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
#' @rdname doEntrez2UniprotMapping
#' @export 
doEntrez2UniprotMapping <- function(entrez.vec){
  # db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep="");
  # db.map <-  readRDS(db.path);
  db.map <-  queryGeneDB("entrez_uniprot", data.org);
  hit.inx <- match(entrez.vec, db.map[, "gene_id"]);
  unips <- db.map[hit.inx, "accession"];    
  mode(unips) <- "character";
  return(unips);
}
