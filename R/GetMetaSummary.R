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
#' @rdname GetMetaSummary
#' @export 
GetMetaSummary<- function(){
  inmex.meta <- readRDS("inmex_meta.rds");
  sel.nms <- unique(inmex.meta$data.lbl)
  sel.nms <- paste(sel.nms, collapse="; ")
  cls.lbls <- unique(inmex.meta$cls.lbl)
  cls.lbls <- paste(cls.lbls, collapse="; ")
  return(c(length(colnames(inmex.meta$data)),nrow(inmex.meta$data), sel.nms, cls.lbls))
}
