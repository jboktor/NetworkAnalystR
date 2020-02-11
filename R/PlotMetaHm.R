#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cmpdNm PARAM_DESCRIPTION
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
#' @rdname PlotMetaHm
#' @export 
#' @importFrom RJSONIO toJSON
PlotMetaHm <-function(cmpdNm){
  
  allmat = readRDS("allmat.rds")
  current.geneset <- readRDS("current_geneset.rds")
  ids <- current.geneset[[cmpdNm]];
  subset = allmat[which(rownames(allmat) %in% ids),]
  if(length(subset)<1){
    subset = allmat[which(rownames(allmat) %in% ids),]
  }
  
  if(inmex.method %in% c("effectsize", "merge")){
    dims = dim(subset)
    rnms = rownames(subset)
    cnms = colnames(subset) 
    m <- mapply(subset, FUN=as.numeric)
    subset <- matrix(data=m, ncol=dims[2], nrow=dims[1])
    rownames(subset) = rnms
    colnames(subset) = cnms
    subset = subset[complete.cases(subset), ];
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
