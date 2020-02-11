#For internal GSEA of meta-analysis of gene sets
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @param netNm PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION
#' @param lib PARAM_DESCRIPTION
#' @param mType PARAM_DESCRIPTION
#' @param BHth PARAM_DESCRIPTION, Default: 0.05
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformMetaPathCombine
#' @export 
PerformMetaPathCombine <- function(name, netNm, method, lib, mType, BHth=0.05){
  if(!performedDE){
    PerformMetaDeAnal();
  }
  library(fgsea)
  curr.geneset <- LoadEnrLib(lib);
  inmex.method <<- "effectsize";
  meta.stat = "null";
  
  allMeta.mat <- readRDS("allMeta.mat.rds");
  #allMeta.mat[allMeta.mat[,2]==0] <- 1e-20
  rankedVec = as.vector(allMeta.mat[,1])*sign(allMeta.mat[,1]);
  
  names(rankedVec) = rownames(allMeta.mat);
  rankedVec = sort(rankedVec)
  rankedVec = rankedVec[unique(names(rankedVec))]
  fgseaRes <- fgsea(pathways = curr.geneset, 
                    stats = rankedVec,
                    minSize=1,
                    maxSize=1000,
                    nperm=5000)
  fgseaRes = fgseaRes[!duplicated(fgseaRes$pathway),]
  fgseaRes = fgseaRes[,c("size","ES", "NES","padj", "pathway", "pval")]
  fgseaRes=fgseaRes[order(-abs(fgseaRes$ES)),]
  fgseaRes=fgseaRes[order(fgseaRes$pval),] 
  
  fgseaRe = data.frame(fgseaRes)
  rownames(fgseaRes) = fgseaRes$pathway
  es.mat = as.matrix(fgseaRes[,c("ES","padj")]);
  colnames(es.mat) = c("EnrichmentScore","Pvalue")
  rownames(es.mat) = fgseaRes$pathway
  sig.inx <- which(es.mat[, "Pvalue"]<=BHth);
  if(length(sig.inx)<10){
    sig.inx = c(1:10)
  }
  metaset.mat <<- es.mat[sig.inx,];
  metaset.mat.all <<- fgseaRes[sig.inx,]
  ii = SetupMetaGSEAStats(name, netNm, BHth, mType, curr.geneset,lib);
  if(ii == 1){
    return(length(sig.inx));
  }else{
    return(0);
  }
}
