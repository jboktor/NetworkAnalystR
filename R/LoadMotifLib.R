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
#' @rdname LoadMotifLib
#' @export 
LoadMotifLib<-function(){
  
  motif.path <- paste(lib.path, data.org, "/motif_set.rds", sep="");
  motif_set<-readRDS(motif.path);
  current.geneset <- motif_set$set;
  set.ids<- names(current.geneset); 
  names(set.ids) <- names(current.geneset) <- motif_set$term;
  current.setlink <<- motif_set$link;
  current.setids <<- set.ids;
  saveRDS(current.geneset, "current_geneset.rds");
  return(current.geneset);
}
