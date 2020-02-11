# perform DE analysis on individual data (w.r.t common matrix)
# to be used/compared in the later analysis, with p-val Inf so that
# de can be adjusted based on user specified in meta later
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param is.meta PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformEachDEAnal
#' @export 
PerformEachDEAnal <- function(is.meta=F){
  inmex.ind <- list();
  inmex.meta <- readRDS("inmex_meta.rds");
  sel.nms <- names(mdata.all)[mdata.all==1];
  if(is.meta){
    for(i in 1:length(sel.nms)){
      dataSet= list()
      dataName <- sel.nms[i];
      sel.inx <- inmex.meta$data.lbl == dataName;
      group <- factor(inmex.meta$cls.lbl[sel.inx]); # note regenerate factor to drop levels 
      data <- inmex.meta$data[, sel.inx];
      
      # update data set
      dataSet$type <- "array";
      dataSet$name <- dataName;
      group <- factor(inmex.meta$cls.lbl[sel.inx]); # note regenerate factor to drop levels 
      dataSet$cls <- group;
      res.limma <- PerformLimma(data, group);
      
      # save the limma fit object for meta-analysis (such as "dataSet1.fit.obj")
      saveRDS(res.limma$fit.obj, file=paste(dataName, "fit.obj", sep="."));
      
      res.all <- GetLimmaResTable(res.limma$fit.obj);
      saveRDS(res.all, "meta.resTable.rds");
      res.mat <- cbind(logFC=res.all$logFC, Pval = res.all$adj.P.Val);
      
      rownames(res.mat) <- rownames(res.all);
      inmex.ind[[dataName]] <- res.mat;
      
      #register sig one
      sig.inx <- res.mat[,2]<=GlobalCutOff$BHth;
      #sig.inx <- res.mat[,2]<=dataSet$pval;
      dataSet$sig.mat <- res.mat[sig.inx,];
      RegisterData(dataSet);
      
      # clean up
      rm(dataSet, res.all);
      gc();
    }
  }else{
    for(i in 1:length(sel.nms)){
      dataName <- sel.nms[i];
      sel.inx <- inmex.meta$data.lbl == dataName;
      group <- factor(inmex.meta$cls.lbl[sel.inx]); # note regenerate factor to drop levels 
      data <- inmex.meta$data[, sel.inx];
      
      dataSet <- readRDS(dataName);
      grp.lvl <- levels(dataSet$cls);
      
      # update data set
      dataSet$type <- "array";
      group <- factor(inmex.meta$cls.lbl[sel.inx], levels=grp.lvl, ordered=T); # note regenerate factor to drop levels 
      dataSet$cls <- group;
      
      res.limma <- PerformLimma(data, group);
      
      # save dataSet object for meta-analysis
      saveRDS(res.limma$fit.obj, file=paste(dataName, "fit.obj", sep="."));
      
      res.all <- GetLimmaResTable(res.limma$fit.obj);
      saveRDS(res.all, "meta.resTable.rds");
      
      res.mat <- cbind(logFC=res.all$logFC, Pval = res.all$adj.P.Val);
      
      rownames(res.mat) <- rownames(res.all);
      inmex.ind[[dataName]] <- res.mat;
      
      #sig.inx <- res.mat[,2]<=GlobalCutOff$BHth;
      sig.inx <- res.mat[,2]<=dataSet$pval;
      dataSet$sig.mat <- res.mat[sig.inx,];
      RegisterData(dataSet);
      # clean up
      rm(dataSet, res.all);
      gc();
    }
  }
  inmex.ind <<- inmex.ind;
}
