#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param my.grps PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetColorSchema
#' @export 
GetColorSchema <- function(my.grps){
  # test if total group number is over 9
  my.grps = as.factor(my.grps);
  grp.num <- length(levels(my.grps));
  
  if(grp.num > 9){
    pal12 = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
              "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
              "#FFFF99", "#B15928");
    dist.cols <- colorRampPalette(pal12)(grp.num);
    lvs <- levels(my.grps);
    colors <- vector(mode="character", length=length(my.grps));
    for(i in 1:length(lvs)){
      colors[my.grps == lvs[i]] <- dist.cols[i];
    }
  }else{
    colors <- as.numeric(my.grps)+1;
  }
  return (colors);
}
