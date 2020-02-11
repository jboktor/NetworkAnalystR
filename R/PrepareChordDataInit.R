#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PrepareChordDataInit
#' @export 
PrepareChordDataInit<-function(){
  # create a list store all possible combination (for a max of 4)
  # note, the whole region are divided into small areas (16 for 4; 7 for 3, 3 for 2)
  # for instance:
  # a: a unique (no b, no c)
  # ab: a and b, no c
  newDat <- list();
  if(anal.type == "metadata"){
    sel.inx <- mdata.all==1;
    sel.nms <- names(mdata.all)[sel.inx];
  }else{
    sel.nms <- listNms;
  }
  venn.genenb <- vector()
  sel.dats <- list();
  for(i in 1:length(sel.nms)){
    nm = sel.nms[i]
    dataSet <- readRDS(nm);
    if(anal.type == "metadata"){
      sel.dats[[nm]] <- rownames(dataSet$sig.mat)
    }else{
      sel.dats[[nm]] <- rownames(dataSet$prot.mat)
    }
    venn.genenb[i] = length(sel.dats[[nm]])
  }
  if(anal.type == "metadata" & meta.selected){
    sel.dats[["meta_dat"]] <- as.character(meta.stat$de);
    venn.genenb[length(venn.genenb) + 1] = length(as.character(meta.stat$de))
  }
  
  chord.list <<- sel.dats;
  chord.genenb <<- venn.genenb
  chord.list.up <<- sel.dats;
  chord.genenb.up <<- venn.genenb
  return(PrepareChordData());
}
