##################################################
## R script for NetworkAnalyst
## Description: prepare data for Venn diagram
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################
# create a list store all possible combination (for a max of 4)
# note, the whole region are divided into small areas (16 for 4; 7 for 3, 3 for 2)
# for instance:
# a: a unique (no b, no c)
# ab: a and b, no c
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
#' @rdname PrepareVennData
#' @export 
PrepareVennData<-function(){
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
  if(length(sel.dats) == 2){
    venn.list <<- Prepare2Venn(sel.dats);
  }else if(length(sel.dats) == 3){
    venn.list <<- Prepare3Venn(sel.dats);
  }else if(length(sel.dats) == 4){
    venn.list <<- Prepare4Venn(sel.dats);
  }else{
    venn.list <<- Prepare4Venn(sel.dats[c(1:4)]);
  }
  venn.list <<- sel.dats;
  venn.genenb <<- venn.genenb
  venn.list.up <<- sel.dats;
  venn.genenb.up <<- venn.genenb
  return(1);
}
