# mapping between genebank, refseq and entrez
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param q.vec PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname doGeneIDMapping
#' @export 
doGeneIDMapping <- function(q.vec, type){
  if(is.null(q.vec)){
    db.map <-  queryGeneDB("entrez", data.org);
    q.vec <- db.map[, "gene_id"];
    type = "entrez";
  }
  if(type == "symbol"){
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "symbol"]);
  }else if(type == "entrez"){
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
  }else{
    # note, some ID can have version number which is not in the database
    # need to strip it off NM_001402.5 => NM_001402
    q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
    q.vec <- q.mat[,1];
    if(type == "genbank"){
      db.map <-  queryGeneDB("entrez_gb", data.org);
    }else if(type == "refseq"){
      db.map <-  queryGeneDB("entrez_refseq", data.org);
    }else if(type == "emblgene"){
      db.map <-  queryGeneDB("entrez_embl_gene", data.org);
    }else if(type == "embltranscript"){
      db.map <-  queryGeneDB("entrez_embl_transcript", data.org);
    }else if(type == "emblprotein"){
      db.map <-  queryGeneDB("entrez_embl_protein", data.org);
    }else if(type == "orfid"){ # only for yeast
      db.map <-  queryGeneDB("entrez_orf", data.org);
      db.path <- paste(lib.path, data.org, "/entrez_orf.rds", sep="");
    }else if(type == "tair"){ # only for ath
      db.map <-  queryGeneDB("tair", data.org);
    }else if(type == "wormbase"){ # only for cel
      db.map <-  queryGeneDB("entrez_wormbase", data.org);
    }else{
      print("Unknown data type");
      return(0);
    }
    hit.inx <- match(q.vec, db.map[, "accession"]);
  }
  entrezs=db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  rm(db.map, q.vec); gc();
  return(entrezs);
}
