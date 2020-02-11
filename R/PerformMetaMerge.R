# This approach directly merge all data sets
# and analyze it as a single data
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param BHth PARAM_DESCRIPTION, Default: 0.05
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformMetaMerge
#' @export 
PerformMetaMerge<-function(BHth=0.05){
  inmex.method <<- "merge";
  meta.mat <<- meta.stat <<- NULL;
  inmex.meta <- readRDS("inmex_meta.rds");
  # prepare for meta-stats
  # calculate sig genes for individual analysis
  shared.names <- rownames(inmex.meta$data);
  
  res.limma <- PerformLimma(inmex.meta$data, as.factor(inmex.meta$cls.lbl));
  res.all <- GetLimmaResTable(res.limma$fit.obj);
  
  ord.inx <- order(res.all$adj.P.Val, decreasing=F);
  dm.mat <- as.matrix(res.all[ord.inx,c("logFC", "adj.P.Val")]);
  colnames(dm.mat) <- c("CombinedLogFC", "Pval");
  
  sig.inx <- which(dm.mat[,"Pval"] <= BHth);
  meta.mat <<- dm.mat[sig.inx,];
  meta.mat.all <<- dm.mat
  SetupMetaStats(BHth);
  return(length(sig.inx));
}
