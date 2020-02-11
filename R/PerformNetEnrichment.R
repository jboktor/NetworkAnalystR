# note: hit.query, resTable must synchronize
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param file.nm PARAM_DESCRIPTION
#' @param fun.type PARAM_DESCRIPTION
#' @param IDs PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformNetEnrichment
#' @export 
PerformNetEnrichment <- function(file.nm, fun.type, IDs){
  # prepare query
  ora.vec <- NULL;
  if(ppi.net$db.type == 'ppi'){
    if(data.org == "ath"){
      idtype <- "tair"
    }else if(data.org %in% c("bsu", "tbr", "cel", "dme", "eco", "pfa", "pae") & net.type == "string"){
      idtype <- "string"
    }else if(data.org %in% c("bta","dre","rno","gga","hsa","mmu") & net.type == "string"){
      idtype <- "emblprotein"
    }else if(data.org %in% c("hsa","mmu", "cel", "dme","sce") & net.type %in% c("innate", "irefinx", "rolland")){
      idtype <- "uniprot"
    }else if(data.org == "sce" & net.type == "string"){ # only for yeast
      idtype <- "emblgene";
    }
    if(idtype=="uniprot"){
      uniprot.vec <- unlist(strsplit(IDs, "; "));
      ora.vec <- doUniprot2EntrezMapping(uniprot.vec);
      names(ora.vec) <- uniprot.vec;
    }else if(idtype=="emblprotein"){
      emblprotein.vec <- unlist(strsplit(IDs, "; "))
      ora.vec <- doEmblProtein2EntrezMapping(emblprotein.vec);
      names(ora.vec) <- emblprotein.vec;
    }else if(idtype=="string"){
      string.vec <- unlist(strsplit(IDs, "; "))
      ora.vec <- doString2EntrezMapping(string.vec);
      names(ora.vec) <- string.vec;
    }else if(idtype=="emblgene"){
      emblgene.vec <- unlist(strsplit(IDs, "; "))
      ora.vec <- doEmblGene2EntrezMapping(emblgene.vec);
      names(ora.vec) <- emblgene.vec;
    }else{
      ora.vec <- unlist(strsplit(IDs, "; "));
      names(ora.vec) <- ora.vec;
    }
  }else{ # net is tf/mir/drug, they already in entrez
    ora.vec <- unlist(strsplit(IDs, "; "));
    names(ora.vec) <- as.character(ora.vec);
  }
  res <- PerformEnrichAnalysis(file.nm, fun.type, ora.vec);
  return(res);
  
}
