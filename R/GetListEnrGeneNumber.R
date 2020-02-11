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
#' @rdname GetListEnrGeneNumber
#' @export 
GetListEnrGeneNumber <- function(){
  all.enIDs <- NULL;
  listSizes <- list();
  if(anal.type == "genelist"){
    if(numOfLists > 1){
      newDat <- list();
      tot.count <- 0;
      all.nms <- listNms;
      for(i in 1:length(all.nms)){
        dataNm <- all.nms[i];
        dataSet <- readRDS(dataNm);
        gene.mat <- dataSet$prot.mat;
        
        # convert to entrez
        expr.val <- gene.mat[,1];
        en.ids <- rownames(gene.mat);
        
        names(expr.val) <- en.ids;
        newDat[[dataNm]] <- expr.val;
        names(en.ids) <- doEntrez2SymbolMapping(en.ids)
        all.enIDs <- c(all.enIDs, en.ids);
        listSizes[[i]] <- list(
          name = dataNm,
          label = dataNm,
          size = length(en.ids)
          #val = de.prct[i]
        )
      }
      
    }else{
      
      all.enIDs <- rownames(dataSet$prot.mat);
      names(all.enIDs ) <- doEntrez2SymbolMapping(all.enIDs)
      listSizes[[1]] = list(
        name = "datalist1",
        label = "datalist1",
        size = length(all.enIDs)
        #val = de.prct[i]
      )
    }
  }else if(anal.type == "onedata"){
    all.enIDs <- rownames(dataSet$sig.mat);
    names(all.enIDs) <- doEntrez2SymbolMapping(all.enIDs)
    listSizes[[1]] = list(
      name = "dataSet1",
      label = "dataSet1",
      size = length(all.enIDs)
      #val = de.prct[i]
    )
  }else{
    newDat <- list();
    tot.count <- 0;
    listSizes <- list();
    all.nms <- names(mdata.all);
    for(i in 1:length(all.nms)){
      dataNm <- all.nms[i];
      dataSet <- readRDS(dataNm);
      gene.mat <- dataSet$sig.mat;
      
      # convert to entrez
      expr.val <- gene.mat[,1];
      en.ids <- rownames(gene.mat);
      
      names(expr.val) <- en.ids;
      newDat[[dataNm]] <- expr.val;
      names(en.ids) <- doEntrez2SymbolMapping(en.ids)
      all.enIDs <- c(all.enIDs, en.ids);
      listSizes[[i]] <- list(
        name = dataNm,
        label = dataNm,
        size = length(en.ids)
      )
    }
  }
  list.genes <<- all.enIDs;
  listSizes <<- listSizes;
}
