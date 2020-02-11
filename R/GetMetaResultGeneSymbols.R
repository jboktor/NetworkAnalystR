# note, due to limitation of get/post
# maximum gene symb for list is top 500
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
#' @rdname GetMetaResultGeneSymbols
#' @export 
GetMetaResultGeneSymbols<-function(){
  ids <- rownames(as.matrix(meta.mat));
  if(length(ids) > 500){
    ids <- ids[1:500];
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  if(inmex.meta$id.type == "entrez"){ # row name gene symbols
    ids <- inmex.meta$gene.symbls[ids];
  }
  return(ids);
}
