#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param gene.id PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PlotSelectedGeneLoading
#' @export 
PlotSelectedGeneLoading<-function(gene.id){
  if(anal.type == "metadata"){
    PlotSelectedGeneMeta(gene.id);
  }else{
    PlotSelectedGene(gene.id, "notVolcano");
  }
}
