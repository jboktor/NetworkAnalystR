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
#' @rdname PrepareSelChordData
#' @export 
PrepareSelChordData<-function(selectedNms){
  newDat <- list();
  sel.nms <- unlist(strsplit(selectedNms, ";"));
  nm.vec <<- sel.nms;
  SelectData();
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
  
  chord.list.up <<- sel.dats;
  chord.genenb.up <<- venn.genenb
  return(PrepareChordData());
}
