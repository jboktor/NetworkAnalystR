# geneIDs is text one string, need to make to vector
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param listNm PARAM_DESCRIPTION
#' @param geneIDs PARAM_DESCRIPTION
#' @param org PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname performGene2ProteinMapping
#' @export 
performGene2ProteinMapping <- function(listNm, geneIDs, org, type){
  
  dataSet <- list();
  dataSet$orig <- geneIDs;
  current.msg <<- NULL;
  data.org <<- org;
  SetInitLib(org)
  
  listNms = vector();
  dataList <- .parseListInput(geneIDs);
  all.prot.mat <- list(); 
  for(i in 1:length(dataList)){
    dataSet$name = paste0("datalist", i);
    listNms[i] = dataSet$name;
    gene.mat <- prot.mat <- dataList[[i]];
    GeneAnotDB <-convertIdToEntrez(rownames(gene.mat), type);
    
    na.inx <- is.na(GeneAnotDB[,1]) | is.na(GeneAnotDB[,2]);
    if(sum(!na.inx) < 2){
      current.msg <<- "Less than two hits found in uniprot database. ";
      print(current.msg);
    }
    rownames(prot.mat) <- GeneAnotDB[,2];
    prot.mat <- prot.mat[!na.inx, , drop=F];
    
    # now merge duplicates
    prot.mat <- RemoveDuplicates(prot.mat, "mean", quiet=T); 
    
    if(is.null(dim(prot.mat))){
      prot.mat <- matrix(prot.mat);
    }
    
    seed.proteins <- rownames(prot.mat);
    dataSet$GeneAnotDB <- GeneAnotDB;
    dataSet$sig.mat <- gene.mat;
    dataSet$prot.mat <- prot.mat;
    dataSet$seeds.proteins <- seed.proteins;
    if(i == 1){
      all.prot.mat <- prot.mat;
      totalseed.proteins = seed.proteins
    }else{
      totalseed.proteins  = c(totalseed.proteins, seed.proteins);
      all.prot.mat <- rbind(all.prot.mat, prot.mat)
    }
    RegisterData(dataSet); 
  }
  all.ent.mat <<- all.prot.mat
  rownames(all.prot.mat) = doEntrez2SymbolMapping(rownames(all.prot.mat))
  all.prot.mat <<- all.prot.mat
  listNms <<- listNms
  mdata.all <- list();
  for(i in 1:length(listNms)){
    mdata.all[i] <- 1;
  }
  names(mdata.all) <- listNms;
  mdata.all <<- mdata.all
  return(totalseed.proteins)
}
