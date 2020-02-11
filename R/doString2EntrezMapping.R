#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param string.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doString2EntrezMapping
#' @export 
doString2EntrezMapping <- function(string.vec){
  # db.path <- paste(lib.path, data.org, "/entrez_string.rds", sep="");
  # db.map <-  readRDS(db.path);
  db.map <-  queryGeneDB("entrez_string", data.org);
  hit.inx <- match(string.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}
