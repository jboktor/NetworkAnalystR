#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param q.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doPpiIDMapping
#' @export 
doPpiIDMapping <- function(q.vec){
  if(data.org == "ath"){
    db.map <-  queryGeneDB("tair", data.org);
  }else if(data.org == "sce"){
    if(net.type == "string"){ # only for yeast
      db.map <-  queryGeneDB("entrez_embl_gene", data.org);
    }else{
      db.map <-  queryGeneDB("entrez_uniprot", data.org);
    }
  }else if(net.type %in% c("innate", "irefinx", "rolland")){
    db.map <-  queryGeneDB("entrez_uniprot", data.org);
  }else{
    if(data.org %in% c("bsu", "tbr", "cel", "dme", "eco", "pfa","pae")){
      db.map <-  queryGeneDB("entrez_string", data.org);
    }else if (data.org %in% c("mmu","hsa") && net.type == "string"){
      db.map <-  queryGeneDB("entrez_string", data.org);
    }else if(data.org == "mtb"){
      db.map <-  queryGeneDB("entrez", data.org);
    }else{
      db.map <-  queryGeneDB("entrez_embl_protein", data.org);
    }
  }
  hit.inx <- match(q.vec, db.map[, "gene_id"]);
  ppi.mat <- db.map[hit.inx, ]; 
  
  # fix the factor col related to library issue
  i <- sapply(ppi.mat, is.factor)
  ppi.mat[i] <- lapply(ppi.mat[i], as.character)
  if(data.org %in% c("pae", "mtb")){
    ppi.mat = ppi.mat[,c(2,1)];
    colnames(ppi.mat) = c("gene_id", "accession");
  }
  return(ppi.mat);
}
