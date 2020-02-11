# Meta-analysis combining effect size
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: 'rem'
#' @param BHth PARAM_DESCRIPTION, Default: 0.05
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformMetaEffectSize
#' @export 
PerformMetaEffectSize<- function(method="rem", BHth=0.05){
  if(!performedDE){
    PerformMetaDeAnal();
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  inmex.method <<- "effectsize";
  meta.mat <<- meta.stat <<- NULL;
  
  sel.nms <- names(mdata.all)[mdata.all==1];
  nbstudies <- length(sel.nms);
  listgd<-vector("list", (nbstudies+3));
  ES<-array(dim=c(nrow(inmex.meta$data),4,nbstudies));
  cls.lvls <- levels(as.factor(inmex.meta$cls.lbl));
  
  for (i in 1:nbstudies){
    data.nm <- sel.nms[i];
    dataSet <- readRDS(data.nm);
    
    fit.obj.nm <- paste(data.nm, "fit.obj", sep=".");
    fit2i <- readRDS(fit.obj.nm);
    
    pvals <- p.adjust(fit2i$p.value,method="BH");
    listgd[[i]]=which(pvals<=BHth);
    
    n1i=length(which(dataSet$cls==cls.lvls[1]));
    n2i=length(which(dataSet$cls==cls.lvls[2]));
    ES[,,i]=effectsize(fit2i$t,((n1i*n2i)/(n1i+n2i)),(fit2i$df.prior+fit2i$df.residual));
  }
  
  #only unbiased; for biased effect sizes, add ES[,1,],ES[,2,]
  listgd[[(nbstudies+1)]]=unique(unlist(listgd[1:nbstudies]))
  restempdirect=combineES(ES[,3,],ES[,4,],BHth, method);
  
  pooled.ef <- restempdirect[[3]];
  wt.mat<- restempdirect[[4]]; # one column for each study
  
  listgd[[(nbstudies+2)]]=restempdirect$DEindices
  listgd[[(nbstudies+3)]]=restempdirect$TestStatistic
  names(listgd)=c(paste("study",1:nbstudies,sep=""),"AllIndStudies","Meta","TestStatistic")  
  
  es.mat <- matrix(0, nrow=nrow(inmex.meta$data), ncol=2);
  es.mat[,1] <- pooled.ef;
  es.mat[,2] <- p.adjust(2*(1-pnorm(abs(listgd[[(nbstudies+3)]]))), method="BH");
  
  colnames(es.mat) <- c("CombinedES","Pval");
  rownames(es.mat) <- rownames(inmex.meta$data);
  #allMeta.mat <<- es.mat;
  saveRDS(es.mat, "allMeta.mat.rds");
  
  # now keep only genes with at least on sig (in one study or meta analysis)
  inx <- union(listgd[[(nbstudies+1)]], listgd[[(nbstudies+2)]]);
  es.mat <- es.mat[inx,];
  
  #sort
  ord.inx <- order(abs(es.mat[, "Pval"]), decreasing = F);
  es.mat<-signif(es.mat[ord.inx,],5);
  
  sig.inx <- which(es.mat[,"Pval"]<=BHth);
  meta.mat <<- es.mat[sig.inx, ];
  meta.mat.all <<- es.mat;
  SetupMetaStats(BHth);
  return(length(sig.inx));
}
