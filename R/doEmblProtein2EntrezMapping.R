#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param emblprotein.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doEmblProtein2EntrezMapping
#' @export 
doEmblProtein2EntrezMapping <- function(emblprotein.vec){
  # db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
  # db.map <-  readRDS(db.path);
  db.map <-  queryGeneDB("entrez_embl_protein", data.org);
  hit.inx <- match(emblprotein.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}
