# perform limma on given two groups selected 
# used by integarative analysis
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataName PARAM_DESCRIPTION
#' @param grps PARAM_DESCRIPTION
#' @param p.lvl PARAM_DESCRIPTION
#' @param fc.lvl PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformLimmaDE
#' @export 
PerformLimmaDE<-function(dataName, grps, p.lvl, fc.lvl=NULL){
  
  print("doing differential analysis ....");
  dataSet <- readRDS(dataName);
  dataSet$pval <- p.lvl
  if(length(levels(dataSet$cls))>2){ 
    grp.nms <- strsplit(grps, " vs. ")[[1]];
    sel.inx <- as.character(dataSet$cls) %in% grp.nms;
  }else{
    sel.inx <- rep(T, ncol(dataSet$data));
  }
  
  group <- factor(dataSet$cls[sel.inx]); # note regenerate factor to drop levels 
  data <- dataSet$data[, sel.inx];
  
  res.limma <- PerformLimma(data, group);
  res.all <- GetLimmaResTable(res.limma$fit.obj);
  
  if(!is.null(fc.lvl)){
    hit.inx <- abs(res.all$logFC)>= fc.lvl & res.all$adj.P.Val <= p.lvl
  }else{
    hit.inx <- res.all$adj.P.Val <= p.lvl
  }
  # note, hit.inx can contain NA, not T/F
  hit.inx <- which(hit.inx);
  res<-res.all[hit.inx,];
  
  # rm .txt suffix for new names
  shortNm <- substring(dataName, 0, nchar(dataName)-4);
  write.csv(signif(res[,-1],5), file=paste("SigGenes_", shortNm, ".csv",sep=""));
  
  sig.count <- nrow(res);
  de.genes <- rownames(res);
  res.mat <- cbind(res.all$logFC, res.all$adj.P.Val);
  rownames(res.mat) <- rownames(res.all);
  non.sig.count <- nrow(data)-sig.count;
  rm(res.all);
  
  gc();
  RegisterData(dataSet);
  # record the sig gene vec
  return (c(1, sig.count, non.sig.count));
}
