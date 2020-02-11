# combining p values based on Fisher's or Stouffer
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: 'stouffer'
#' @param BHth PARAM_DESCRIPTION, Default: 0.05
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformPvalCombination
#' @export 
PerformPvalCombination <- function(method="stouffer", BHth=0.05){
  if(!performedDE){
    PerformMetaDeAnal();
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  inmex.method <<- "metap";
  meta.mat <<- meta.stat <<- NULL;
  sel.nms <- names(mdata.all)[mdata.all==1];
  
  classes <- list();
  nbstudies <- length(sel.nms);
  listgd=vector("list", (nbstudies+3));
  
  for (i in 1:nbstudies){
    data.nm <- sel.nms[i];
    dataSet <- readRDS(data.nm);
    classes[[i]] <- dataSet$cls; 
    
    fit.obj.nm <- paste(data.nm, "fit.obj", sep=".");
    fit2i <- readRDS(fit.obj.nm);
    
    pvals <- p.adjust(fit2i$p.value,method="BH");
    listgd[[i]]=which(pvals<=BHth);
    
    #recalculate moderated one sided p
    p1sidedLimma=pt(fit2i$t,df=(fit2i$df.prior+fit2i$df.residual))
    assign(paste("p1sidedLimma",i,sep=""), p1sidedLimma)
  }
  
  names(classes) <- sel.nms;
  tempvec=paste("p1sidedLimma",1:nbstudies,sep="");
  
  lsinglep=lapply(tempvec,FUN=function(x) get(x,inherits=TRUE));
  nrep=unlist(lapply(classes,FUN=function(x)length(x)));
  listgd[[(nbstudies+1)]]=unique(unlist(listgd[1:nbstudies]));
  
  restempdirect=combinePvals(lsinglep,nrep,BHth,method);
  
  listgd[[(nbstudies+2)]]=restempdirect$DEindices
  listgd[[(nbstudies+3)]]=restempdirect$CombinedP
  names(listgd)=c(paste("study",1:nbstudies,sep=""),"AllIndStudies","Meta","CombinedP");  
  
  pc.mat <- cbind(CombinedTstat=restempdirect$CombinedStat, CombinedPval=restempdirect$CombinedP);
  rownames(pc.mat) <- rownames(inmex.meta$data);
  saveRDS(pc.mat, "allMeta.mat.rds");
  
  # now keep only genes with at least on sig (in one study or meta analysis)
  inx <- union(listgd[[(nbstudies+1)]], listgd[[(nbstudies+2)]]);
  pc.mat <- pc.mat[inx,];
  
  #sort
  ord.inx <- order(pc.mat[, "CombinedPval"], decreasing = F);
  pc.mat<-signif(pc.mat[ord.inx,],5);
  
  sig.inx <- which(pc.mat[, "CombinedPval"]<=BHth);
  meta.mat <<- pc.mat[sig.inx, ];
  meta.mat.all <<- pc.mat;
  SetupMetaStats(BHth);
  return(length(sig.inx));
}
