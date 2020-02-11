#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cmpdNm PARAM_DESCRIPTION
#' @param IDs PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[RJSONIO]{toJSON}}
#' @rdname PlotGShm
#' @export 
#' @importFrom RJSONIO toJSON
PlotGShm <-function(cmpdNm, IDs){
  ids = unlist(strsplit(IDs, "; "));
  cmpdNm <- gsub(" ", "_",  cmpdNm);
  cmpdNm <- gsub("/", "_",  cmpdNm);
  if(anal.type == "onedata"){
    subset = dataSet$data.norm[which(doEntrez2SymbolMapping(rownames(dataSet$data.norm)) %in% ids),]
    if(length(subset)<1){
      subset = dataSet$data.norm[which(rownames(dataSet$data.norm) %in% ids),]
    }
  }else{
    if(dataSet$name != selDataNm){
      dataSet <- readRDS(selDataNm);
    }
    subset = dataSet$data[which(doEntrez2SymbolMapping(rownames(dataSet$data)) %in% ids),]
    if(length(subset)<1){
      subset = dataSet$data[which(rownames(dataSet$data) %in% ids),]
    }
  }
  
  json.res <- list(
    data=subset,
    ids=doEntrez2SymbolMapping(rownames(subset)),
    entrez = rownames(subset)
  )
  library(RJSONIO);
  json.mat <- RJSONIO::toJSON(json.res, .na='null', pretty=TRUE);
  json.nm <- paste(cmpdNm,"_hm", ".json", sep="");
  
  sink(json.nm);
  cat(json.mat);
  sink();
  
  return(json.nm)
}
