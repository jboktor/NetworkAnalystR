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
#' @rdname PlotGShmMeta
#' @export 
#' @importFrom RJSONIO toJSON
PlotGShmMeta <-function(cmpdNm, IDs){
  ids = unlist(strsplit(IDs, "; "));
  cmpdNm <- gsub(" ", "_",  cmpdNm);
  cmpdNm <- gsub("/", "_",  cmpdNm);
  inmex.meta <- readRDS("inmex_meta.rds");
  subset = dataSet$data.norm[which(doEntrez2SymbolMapping(rownames(inmex.meta$data)) %in% ids),]
  if(length(subset)<1){
    subset = dataSet$data.norm[which(rownames(inmex.meta$data) %in% ids),]
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
