#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param listNm PARAM_DESCRIPTION
#' @param org PARAM_DESCRIPTION
#' @param geneIDs PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformListAnnot
#' @export 
PerformListAnnot <- function(listNm, org, geneIDs, type){
  dataSet <- list();
  dataSet$orig <- "";
  current.msg <<- NULL;
  data.org <<- org;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  SetInitLib(org)
  listNms <- multiFileNamesU;
  numOfLists <<-length(multiFileNamesU);
  notOk = 0
  for(i in 1:length(listNms)){
    dataSet = readRDS(listNms[i])
    dataSet$name = listNms[i];
    gene.mat <- prot.mat <- dataSet$prot.mat;
    GeneAnotDB <-convertIdToEntrez(rownames(gene.mat), type);
    
    na.inx <- is.na(GeneAnotDB[,1]) | is.na(GeneAnotDB[,2]);
    if(sum(!na.inx) < 2){
      current.msg <<- paste0("Less than two hits found in database for ", listNms[i]);
      print(current.msg);
      return(0);
    }
    rownames(prot.mat) <- GeneAnotDB[,2];
    prot.mat <- prot.mat[!na.inx, , drop=F];
    
    # now merge duplicates
    listInxU <<- listNms[i];
    prot.mat <- RemoveDuplicates(prot.mat, "mean", quiet=F); 
    
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
