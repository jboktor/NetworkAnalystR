# perform differential analysis
# default: all pair-wise comparison (A-B) + (B-C) + (A-C)
# custom: only compare two groups (A-C)
# time: compare consecutive groups (B-A) + (C-B)
# reference: all others against common reference (A-C) + (B-C)
# nested: (A-B)+(C-D) 
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataName PARAM_DESCRIPTION
#' @param anal.type PARAM_DESCRIPTION
#' @param par1 PARAM_DESCRIPTION, Default: NULL
#' @param par2 PARAM_DESCRIPTION, Default: NULL
#' @param nested.opt PARAM_DESCRIPTION, Default: 'intonly'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformIndLimmaDE
#' @export 
PerformIndLimmaDE<-function(dataName, anal.type, par1=NULL, par2=NULL, nested.opt = "intonly"){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  require(limma);
  current.msg <<- "";
  cls <- dSet$cls; 
  
  cls.lvls <- levels(cls);
  data <- dSet$data.norm;
  meta.info <- dSet$meta.info;
  
  design <- model.matrix(~ -1 + cls) # no intercept
  colnames(design) <- cls.lvls;
  
  if(is.null(dSet$block)){
    fit = lmFit(data, design);
  }else{
    # limma user guide chapter 8 P49 
    corfit<-duplicateCorrelation(data,design,block=dSet$block);
    fit = lmFit(data, design, block=dSet$block, correlation=corfit$consensus);
  }
  # sanity check
  if(!is.fullrank(design)){
    current.msg <<- paste("This metadata combination is not full rank! Please use other combination."); 
    return(0);
  }
  df.residual <- fit$df.residual;
  if (all(df.residual == 0)){
    current.msg <<- paste("There is not enough replicates in each group (no residual degrees of freedom)!"); 
    return(0);
  }
  
  myargs <- list();
  grp.nms <- cls.lvls;
  
  if(anal.type == 'default'){ # all pair-wise
    inx = 0;
    for(m in 1:(length(grp.nms)-1)){
      for(n in (m+1):length(grp.nms)){
        inx <- inx + 1;
        myargs[[inx]] <- paste(grp.nms[m], "-", grp.nms[n], sep="")
      }
    }
  }else if(anal.type == 'time'){
    for(i in 2:length(grp.nms)){
      myargs[[i-1]] <- paste(grp.nms[i], "-", grp.nms[i-1], sep="")
    }
  }else if(anal.type == 'custom'){
    grp.nms <- strsplit(par1, " vs. ")[[1]];
    myargs[[1]] <- paste(grp.nms, collapse="-");
  }else if(anal.type == 'reference'){
    ref <- par1;
    cntr.cls <- grp.nms[grp.nms!= ref];
    myargs <- as.list(paste(cntr.cls, "-", ref, sep=""));
  }else if(anal.type == 'nested'){
    grp.nms1 <- strsplit(par1, " vs. ")[[1]];
    grp.nms2 <- strsplit(par2, " vs. ")[[1]];
    if(all(grp.nms1 == grp.nms2)){
      current.msg <<- paste("The two nested groups are the same. Please choose two different groups."); 
      return(0);
    }
    grp.nms <- unique(c(grp.nms1, grp.nms2));
    if(nested.opt == "intonly"){
      myargs[[1]] <- paste("(", paste(grp.nms1, collapse="-"), ")-(", paste(grp.nms2, collapse="-"), ")", sep=""); 
    }else{
      myargs[[1]] <- paste(grp.nms1, collapse="-");
      myargs[[2]] <- paste(grp.nms2, collapse="-"); 
      myargs[[3]] <- paste("(", paste(grp.nms1, collapse="-"), ")-(", paste(grp.nms2, collapse="-"), ")", sep=""); 
    }
  }else{ # 
    print(paste('Not supported: ', anal.type));
  }
  
  myargs[["levels"]] <- design;
  contrast.matrix <- do.call(makeContrasts, myargs);
  fit2 <- contrasts.fit(fit, contrast.matrix);
  fit2 <- eBayes(fit2);
  
  resTable <- topTable(fit2, number=Inf, adjust.method="fdr");
  dSet$sig.orig <- resTable; # record orignal data for update
  anot.id <- rownames(resTable);
  gene.anot <- doEntrezIDAnot(anot.id);
  current.msg <<- current.msg; 
  
  # may need to update data and class
  hit.inx <- cls %in% grp.nms;
  if(sum(hit.inx) < length(hit.inx)){
    current.msg <<- paste(current.msg, "Only groups selected for comparisons: ", paste(grp.nms, collapse=", "), "are included.");
    cls <- factor(cls[hit.inx]);
    cls.lvls <- levels(cls);
    data <- data[,hit.inx];
    meta.info <- dSet$meta.info[hit.inx,];
  }
  
  dSet$data.stat <- data;
  dSet$cls.stat <- cls;
  dSet$meta.stat <- meta.info;
  dSet$sig.genes.anot <- gene.anot;
  RegisterData(dataSet);
  return (1);
}
