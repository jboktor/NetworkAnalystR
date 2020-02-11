#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param uniprot.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doUniprot2EntrezMapping
#' @export 
doUniprot2EntrezMapping <- function(uniprot.vec){
  # db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep="");
  # db.map <-  readRDS(db.path);
  db.map <-  queryGeneDB("entrez_uniprot", data.org);
  hit.inx <- match(uniprot.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];    
  mode(entrezs) <- "character";
  return(entrezs);
}
