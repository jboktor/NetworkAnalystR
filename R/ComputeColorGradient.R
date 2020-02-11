#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nd.vec PARAM_DESCRIPTION
#' @param background PARAM_DESCRIPTION, Default: 'black'
#' @param centered PARAM_DESCRIPTION
#' @param colorblind PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ComputeColorGradient
#' @export 
ComputeColorGradient <- function(nd.vec, background="black", centered, colorblind){
  require("RColorBrewer");
  
  minval = min(nd.vec, na.rm=TRUE);
  maxval = max(nd.vec, na.rm=TRUE);
  res = maxval-minval;
  
  if(res == 0){
    return(rep("#FF0000", length(nd.vec)));
  }
  color <- GetColorGradient(background, centered, colorblind);
  breaks <- generate_breaks(nd.vec, length(color), center = centered);
  return(scale_vec_colours(nd.vec, col = color, breaks = breaks));
}
