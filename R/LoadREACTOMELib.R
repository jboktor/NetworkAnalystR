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
#' @rdname LoadREACTOMELib
#' @export 
LoadREACTOMELib<-function(){
  
  reactome.path <- paste(lib.path, data.org, "/reactome.rds", sep="");
  reactome.anot <- readRDS(reactome.path)
  current.geneset <- reactome.anot$sets;
  set.ids<- names(current.geneset); 
  names(set.ids) <- names(current.geneset) <- reactome.anot$term;
  current.setlink <<- reactome.anot$link;
  current.setids <<- set.ids;
  saveRDS(current.geneset, "current_geneset.rds");
  return(current.geneset);
}
