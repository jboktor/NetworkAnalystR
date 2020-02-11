#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param selectedNms PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PrepareSelVennData
#' @export 
PrepareSelVennData<-function(selectedNms){
  newDat <- list();
  sel.nms <- unlist(strsplit(selectedNms, ";"));
  venn.genenb <- vector()
  sel.dats <- list();
  for(i in 1:length(sel.nms)){
    nm = sel.nms[i]
    if(nm != "meta_dat"){
      dataSet <- readRDS(nm);
      if(anal.type == "metadata"){
        sel.dats[[nm]] <- rownames(dataSet$sig.mat)
      }else{
        sel.dats[[nm]] <- rownames(dataSet$prot.mat)
      }
      venn.genenb[i] = length(sel.dats[[nm]])
    }else{
      sel.dats[[nm]] <- as.character(meta.stat$de);
      venn.genenb[i] = length(as.character(meta.stat$de))
    }
  }
  if(length(sel.dats) == 2){
    venn.list <- Prepare2Venn(sel.dats);
  }else if(length(sel.dats) == 3){
    venn.list <- Prepare3Venn(sel.dats);
  }else if(length(sel.dats) == 4){
    venn.list <- Prepare4Venn(sel.dats);
  }
  venn.list.up <<- sel.dats;
  venn.genenb.up <<- venn.genenb
  return(1);
}
