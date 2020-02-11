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
#' @rdname PrepareListChordData
#' @export 
PrepareListChordData <- function(){
  
  all.enIDs <- NULL;
  newDat <- list();
  tot.count <- 0;
  all.nms <- names(mdata.all)[mdata.all==1];
  for(i in 1:length(all.nms)){
    dataNm <- all.nms[i];
    dataSet <- readRDS(dataNm);
    gene.mat <- dataSet$prot.mat;
    
    # convert to entrez
    expr.val <- gene.mat[,1];
    en.ids <- rownames(gene.mat);
    
    names(expr.val) <- en.ids;
    newDat[[dataNm]] <- expr.val;
    
    all.enIDs <- c(all.enIDs, en.ids);
    tot.count <- tot.count + nrow(gene.mat);
    
    if(tot.count > 2000){
      current.msg <<- paste("Chord diagrams is effective to display relationships for less than 1000 items. The results contain", tot.count, 
                            "of genes (max. allowed: 2000). You can try Venn diagram instead.")
      return(NULL);
    }
  }
  PrepareChordDataFromList(newDat, unique(all.enIDs));
}
