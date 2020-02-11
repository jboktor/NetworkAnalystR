#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param c1 PARAM_DESCRIPTION, Default: 'grey'
#' @param c2 PARAM_DESCRIPTION, Default: 'red'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname color_scale
#' @export 
color_scale <- function(c1="grey", c2="red") {
  pal <- colorRampPalette(c(c1, c2))
  colors <- pal(100)
  return(colors)
}
