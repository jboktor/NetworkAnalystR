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
#' @rdname GetMetaResultGeneIDLinks
#' @export 
GetMetaResultGeneIDLinks <- function(){
  ids <- rownames(as.matrix(meta.mat));
  if(length(ids) > 500){
    ids <- ids[1:500];
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  symbs <- inmex.meta$gene.symbls[ids];
  # set up links to genbank
  annots <- paste("<a href='http://www.ncbi.nlm.nih.gov/gene?term=", ids,
                  "' target='_blank'>", symbs, "</a>", sep="");
  return(annots);
}
