#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param emblgene.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doEmblGene2EntrezMapping
#' @export 
doEmblGene2EntrezMapping <- function(emblgene.vec){
  # db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
  # db.map <-  readRDS(db.path);
  db.map <-  queryGeneDB("entrez_embl_gene", data.org);
  hit.inx <- match(emblgene.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}
