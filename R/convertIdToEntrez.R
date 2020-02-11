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
#' @rdname convertIdToEntrez
#' @export 
convertIdToEntrez <- function(q.vec, type){ #convert user input ids to entrez
  if(type == "entrez"){
    # need to get only our data
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
    
  }else if(type == "symbol"){
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "symbol"]);
    
  }else{
    if(type == "genbank"){
      # note, some ID can have version number which is not in the database
      # need to strip it off NM_001402.5 => NM_001402
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.map <-  queryGeneDB("entrez_gb", data.org);
    }else if(type == "refseq"){
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.map <-  queryGeneDB("entrez_refseq", data.org);
    }else if(type == "emblgene"){
      db.map <-  queryGeneDB("entrez_embl_gene", data.org);
    }else if(type == "tair"){
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.map <-  queryGeneDB("tair", data.org);
    }else if(type == "embltranscript"){
      db.map <-  queryGeneDB("entrez_embl_transcript", data.org);
    }else if(type == "emblprotein"){
      db.map <-  queryGeneDB("entrez_embl_protein", data.org);
    }else if(type == "orfid"){ # only for yeast
      db.map <-  queryGeneDB("entrez_orf", data.org);
    }else if(type == "flybase"){
      db.map <-  queryGeneDB("entrez_flybase", data.org);
    }else if(type == "string"){ 
      db.map <-  queryGeneDB("entrez_string", data.org);
    }else if(type == "ecogene"){ # only for ecoli
      db.map <-  queryGeneDB("entrez_ecogene", data.org);
    }else if(type == "uniprot"){
      db.map <-  queryGeneDB("entrez_uniprot", data.org);
    }else if(type == "paelocus"){
      db.map <-  queryGeneDB("entrez", data.org);
    }else{
      print("Unknown data type");
      return(0);
    }
    hit.inx <- match(q.vec, db.map[, "accession"]);
    entrezs <- db.map[hit.inx, ];
  }
  entrezs <- db.map[hit.inx, ]; 
  if(type == "entrez"){
    entrezs = entrezs[,c(1,1)];
  }else{
    entrezs = entrezs[,c(2,1)];
  }
  colnames(entrezs) <- c("accession", "gene_id")
  return(entrezs);
}
