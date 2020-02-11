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
#' @rdname GetMetaResultPathLinks
#' @export 
GetMetaResultPathLinks <- function(){
  symbs <- rownames(meta.mat);
  ids <- current.setids[symbs];
  # set up links to genbank
  annots <- paste("<a href='http://pantherdb.org/panther/category.do?categoryAcc=", ids,
                  "' target='_blank'>", symbs, "</a>", sep="");
  return(annots);
}
