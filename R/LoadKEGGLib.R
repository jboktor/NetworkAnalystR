# Load various libaries for functional enrichment analysis
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
#' @rdname LoadKEGGLib
#' @export 
LoadKEGGLib<-function(){
  kegg.path <- paste(lib.path, data.org, "/kegg.rds", sep="");
  
  kegg.anot <- readRDS(kegg.path)
  current.setlink <- kegg.anot$link;
  current.geneset <- kegg.anot$sets;
  set.ids<- names(current.geneset); 
  names(set.ids) <- names(current.geneset) <- kegg.anot$term;
  
  current.setlink <<- current.setlink;
  current.setids <<- set.ids;
  saveRDS(current.geneset, "current_geneset.rds");
  return(current.geneset);
}
