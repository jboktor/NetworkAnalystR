##################################################
## R scripts for NetworkAnalyst 
## Description: Meta Analysis Methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# for multiple class, only select two
# also record all grp lbls
SetGroupContrast <- function(dataName, grps){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  if(length(levels(dataSet$cls))>2){ 
    print("Updating group contrasts .....");
    grp.nms <- strsplit(grps, " vs. ")[[1]];
    sel.inx <- as.character(dataSet$cls) %in% grp.nms;
    
    # regenerate factor to drop levels, force the levels order
    group <- factor(dataSet$cls[sel.inx], levels=grp.nms);  
    data <- dataSet$data[, sel.inx];
    dataSet$cls <- group;
    dataSet$data <- data;
    RegisterData(dataSet);  
  }
}

# determine if all annotated data are ready for meta-analysis
CheckMetaDataIntegrity<-function(){
  
  performedDE <<- FALSE;
  
  if(length(mdata.all) == 0){
    current.msg <<-"Please upload your data or try our example datasets!";
    print(current.msg);
    return(0);
  }
  
  include.inx <- mdata.all==1;
  if(sum(include.inx) == 0){
    current.msg <<-"No dataset is selected for analysis!";
    print(current.msg);
    return(0);
  }
  
  sel.nms <- names(mdata.all)[include.inx];
  clss <- list();
  if(meta.upload){
    # update meta data only for select/deselect datasets
    inmex.meta.orig <- readRDS("inmex.meta.orig.rds");
    hit.inx <- inmex.meta.orig$data.lbl %in% sel.nms;
    data <- inmex.meta.orig$data[, hit.inx];
    id.type <- inmex.meta.orig$id.type;
    cls.lbl <- factor(inmex.meta.orig$cls.lbl[hit.inx]);
    data.lbl <- factor(inmex.meta.orig$data.lbl[hit.inx]);
    common.matrix <- data;
    nms.vec <<- rownames(inmex.meta.orig$data);
    smps.vec <<- data.lbl
  }else{   
    # first check that all class labels are consistent
    dataSet <- readRDS(sel.nms[1]);
    lvls <- levels(dataSet$cls);
    id.type <- dataSet$id.type;
    clss[[1]] = dataSet$cls;
    nms <- rownames(dataSet$data);
    shared.nms <- nms;
    for(i in 2:length(sel.nms)){
      dataSet <- readRDS(sel.nms[i]);
      clss[[i]] = dataSet$cls;
      # check if class label is consistent
      if(!all(levels(dataSet$cls) == lvls)){
        current.msg <<- paste(sel.nms[i], "has different group labels", 
                              paste(levels(dataSet$cls), collapse=":"), 
                              "from", sel.nms[1], paste(lvls, collapse=":"));
        print(current.msg);
        return(0);
      }
      
      # check and record if there is common genes            
      shared.nms <- intersect(shared.nms, rownames(dataSet$data));
      if(length(shared.nms) < 10){
        current.msg <<- paste(sel.nms[i], "has less than 10 common genes/probes from previous data sets");
        print(current.msg);
        return(0);
      }
      
      # check gene id type
      if(dataSet$id.type != id.type){
        current.msg <<- paste(sel.nms[i], "has different gene/probe ID from", sel.nms[1]);
        print(current.msg);
        return(0);
      }
    }
    
    nrepu<<-unlist(lapply(clss,FUN=function(x)length(x)));      
    
    print("Passed exp condition check!");
    
    # now construct a common matrix to faciliated plotting across all studies
    dataName <- sel.nms[1];
    dataSet <- readRDS(dataName);
    common.matrix <- dataSet$data[shared.nms, ];
    nms.vec = rownames(dataSet$data);
    smps.vec = colnames(dataSet$data);
    data.lbl <- rep(dataName, ncol(common.matrix));
    cls.lbl <- dataSet$cls;
    
    for(i in 2:length(sel.nms)){
      dataName <- sel.nms[i];
      dataSet <- readRDS(dataName);
      ndat <- dataSet$data[shared.nms, ];
      nms.vec = c(nms.vec, rownames(dataSet$data));
      smps.vec = c(smps.vec, colnames(dataSet$data));
      plot.ndat <- t(scale(t(ndat)));
      common.matrix <- cbind(common.matrix, ndat);
      data.lbl <- c(data.lbl, rep(dataName, ncol(dataSet$data[,])));
      cls.lbl <- c(cls.lbl, dataSet$cls);
    }
    cls.lbl <- factor(cls.lbl);
    levels(cls.lbl) <- lvls;
    
    smps.nms = colnames(common.matrix)
    if(length(unique(smps.nms)) != length(smps.nms)){
      data.nb = length(unique(data.lbl));
      data.vec = vector()
      for(i in 1:data.nb){
        data.vec[i] = paste0("d", i);
      }
      levels(data.lbl) = data.vec;
      colnames(common.matrix) = make.unique(paste(data.vec, smps.nms, sep="_"));
      
      dataSet <- readRDS(sel.nms[1]);
      colnames(dataSet$data) = paste("d1", colnames(dataSet$data), sep="_");
      RegisterData(dataSet);
      
      for(i in 2:length(sel.nms)){
        dataSet <- readRDS(sel.nms[i]);
        colnames(dataSet$data) = paste0("d",i,"_",colnames(dataSet$data));
        # check if class label is consistent
        RegisterData(dataSet);
      }
      smps.vec <<- smps.nms;
      current.msg <<- paste("Duplicated sample names detected, samples have been renamed to make them unique.");
    }
    
    # note: index by entrez, gene symbol DO have duplicates
    rownames(common.matrix) <- shared.nms;
    
    # resort data, first on data.lbl, then on class lbl
    ord.inx <- order(data.lbl, cls.lbl);
    common.matrix <- data.matrix(common.matrix[,ord.inx]);
    cls.lbl <- cls.lbl[ord.inx];
    data.lbl <- data.lbl[ord.inx];
    smps.vec <- smps.vec[ord.inx];
    nms.vec = unique(nms.vec)
    nms.vec <<- nms.vec
    smps.vec <<- smps.vec
  }
  
  if(ncol(common.matrix) > 1000){  # max sample number allow 1000
    current.msg <<- paste("Total combined sample #:", ncol(common.matrix), "(exceed the limit: 1000!)");
    return(0);
  }
  
  # save the meta-dataset
  res <- data.frame(colnames(common.matrix), cls.lbl, data.lbl, t(common.matrix));
  colnames(res) <- c('#NAME', '#CLASS.condition', paste('#CLASS.dataset',id.type, data.org, sep="."), rownames(common.matrix));
  write.table(t(res), sep="\t", file="NetworkAnalyst_merged_data.txt", col.names=F, quote=FALSE);
  
  # need to set up the data for plotting (boxplot, heatmpa) so 
  # we need to scale row for each dataset in order to elimiate the maganitude difference 
  plot.matrix <- matrix(NA, nrow=nrow(common.matrix), ncol=ncol(common.matrix));
  rownames(plot.matrix) <- rownames(common.matrix);
  colnames(plot.matrix) <- colnames(common.matrix);
  for(i in 1:length(sel.nms)){
    data.inx <- data.lbl == sel.nms[i];
    plot.matrix[,data.inx] <- t(scale(t(common.matrix[,data.inx])));
  }
  
  # if entrez, get symbols for display
  shared.nms <- rownames(common.matrix);
  if(id.type == "entrez"){ 
    symbols <- doEntrez2SymbolMapping(shared.nms);
  }else{ # display itself
    symbols <- shared.nms;
  }
  names(symbols) <- shared.nms;
  
  inmex.meta <- list(data=common.matrix,
                     plot.data=plot.matrix,
                     id.type = id.type,
                     gene.symbls = symbols,
                     cls.lbl=factor(cls.lbl),
                     data.lbl=data.lbl);
  
  saveRDS(inmex.meta, "inmex_meta.rds");
  smps.vec <<- colnames(common.matrix);
  
  # setup common stats gene number, smpl number, grp info
  current.msg <<- paste("Sample #:", ncol(inmex.meta$data),
                        "Common ID #:", nrow(inmex.meta$data), 
                        "Condition:", paste(levels(inmex.meta$cls.lbl), collapse=" vs. "));
  return(1);
}

PerformMetaDeAnal <- function(){ 
  inmex.meta <- readRDS("inmex_meta.rds");   
  data.lbl = inmex.meta$data.lbl
  allmat = matrix("NA", nrow=length(nms.vec), ncol=length(data.lbl))
  rownames(allmat) = nms.vec
  colnames(allmat) = colnames(inmex.meta$data)
  
  include.inx <- mdata.all==1;
  if(sum(include.inx) == 0){
    current.msg <<-"No dataset is selected for analysis!";
    print(current.msg);
    return(0);
  }
  
  if(meta.upload){
    dat = inmex.meta$data 
    datasets = levels(inmex.meta$data.lbl)
    for(i in 1:length(datasets)){
      dt = datasets[i];
      inx = which(inmex.meta$data.lbl == dt)
      ind_dat = inmex.meta$data[,inx]
      cols <- colnames(allmat)[colnames(allmat) %in% colnames(ind_dat)]
      rows <- rownames(allmat)[rownames(allmat) %in% rownames(ind_dat)]
      norm.dat = t(apply(ind_dat, 1, function(x){as.numeric(cut(x, breaks=30))}))
      norm.dat = as.matrix(norm.dat)
      colnames(norm.dat) = colnames(ind_dat)
      allmat[rows, cols] <- norm.dat[rows, cols]
    }
  }else{
    sel.nms <- names(mdata.all)[include.inx];
    for(i in 1:length(sel.nms)){
      dataName <- sel.nms[i];
      dataSet <- readRDS(dataName);
      dataSet$data.orig <- NULL;
      cols <- colnames(allmat)[colnames(allmat) %in% colnames(dataSet$data)]
      rows <- rownames(allmat)[rownames(allmat) %in% rownames(dataSet$data)]
      norm.dat = t(apply(dataSet$data, 1, function(x){as.numeric(cut(x, breaks=30))}))
      norm.dat = as.matrix(norm.dat)
      colnames(norm.dat) = colnames(dataSet$data)
      allmat[rows, cols] <- norm.dat[rows, cols]
    }
  }
  
  saveRDS(allmat, "allmat.rds");
  performedDE <<- TRUE;
  PerformEachDEAnal(meta.upload);
}


# perform DE analysis on individual data (w.r.t common matrix)
# to be used/compared in the later analysis, with p-val Inf so that
# de can be adjusted based on user specified in meta later
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

# Meta-analysis combining effect size
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

# combining p values based on Fisher's or Stouffer
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

# diff used for direction, not selection
PerformVoteCounting <- function(BHth = 0.05, minVote){
  if(!performedDE){
    PerformMetaDeAnal();
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  inmex.method <<- "votecount";
  DE.vec <<- NULL; # store entrez id from meta-analysis for GO
  meta.mat <<- meta.stat <<- NULL;
  sel.nms <- names(mdata.all)[mdata.all==1];
  # first create a matrix to stall the result
  # row for each feature and col for each dataset uploaded
  vc.mat <- matrix(0, nrow=nrow(inmex.meta$data), ncol=length(sel.nms)+1);
  shared.ids <- rownames(inmex.meta$data);
  for(i in 1:length(inmex.ind)){
    res.mat <- inmex.ind[[i]];
    res.mat <- res.mat[shared.ids, ];
    
    #note in meta-analysis should consider directions
    # use logFC for this purpose 
    # consider upregulated
    hit.up.inx <- res.mat[,1]> 0 & res.mat[,2] <= BHth;
    up.vote <- as.numeric(hit.up.inx);
    
    # consider downregulated
    hit.dn.inx <- res.mat[,1] < 0 & res.mat[,2] <= BHth;
    dn.vote <- -as.numeric(hit.dn.inx);
    
    vc.mat[,i] <- up.vote + dn.vote;
  }
  
  # total score (votes for each direction)
  vc.mat[,length(sel.nms)+1] <- apply(vc.mat, 1, sum);
  colnames(vc.mat) <- c(paste("Vote", substring(sel.nms,0, nchar(sel.nms)-4)), "VoteCounts");
  rownames(vc.mat) <- rownames(inmex.meta$data);
  
  # compute at least one vote (no direction)
  vote.any <- apply(abs(vc.mat), 1, sum)
  vote.any.inx <- vote.any > 0;
  
  # return results with at least one vote
  vc.mat <- vc.mat[vote.any.inx, ];
  
  #sort
  ord.inx <- order(abs(vc.mat[, "VoteCounts"]), decreasing = T);
  vc.mat <- vc.mat[ord.inx, "VoteCounts", drop=F];
  
  sig.inx <- abs(vc.mat[,"VoteCounts"]) >= minVote;
  meta.mat <<- vc.mat;
  meta.mat.all <<- vc.mat;
  SetupMetaStats(BHth);
  return(length(sig.inx));
}


# This approach directly merge all data sets
# and analyze it as a single data
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

GetMetaGeneIDType<-function(){
  inmex.meta <- readRDS("inmex_meta.rds");
  return(inmex.meta$id.type);
}

GetMetaResultGeneIDs<-function(){
  rnms <- rownames(as.matrix(meta.mat));# already sorted based on meta-p values
  if(length(rnms) > 500){
    rnms <- rnms[1:500];
  }
  return(rnms);
}

# note, due to limitation of get/post
# maximum gene symb for list is top 500
GetMetaResultGeneSymbols<-function(){
  ids <- rownames(as.matrix(meta.mat));
  if(length(ids) > 500){
    ids <- ids[1:500];
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  if(inmex.meta$id.type == "entrez"){ # row name gene symbols
    ids <- inmex.meta$gene.symbls[ids];
  }
  return(ids);
}

GetMetaResultPathSymbols<-function(){
  return(rownames(meta.mat.all));
}

GetMetaResultPathLinks <- function(){
  symbs <- rownames(meta.mat);
  ids <- current.setids[symbs];
  # set up links to genbank
  annots <- paste("<a href='http://pantherdb.org/panther/category.do?categoryAcc=", ids,
                  "' target='_blank'>", symbs, "</a>", sep="");
  return(annots);
}


GetMetaResultGeneIDLinks <- function(){
  ids <- rownames(as.matrix(meta.mat));
  if(length(ids) > 500){
    ids <- ids[1:500];
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  symbs <- inmex.meta$gene.symbls[ids];
  # set up links to genbank
  annots <- paste("<a href='http://www.ncbi.nlm.nih.gov/gene?term=", ids,
                  "' target='_blank'>", symbs, "</a>", sep="");
  return(annots);
}

GetMetaResultColNames<-function(){
  sel.nms <- names(mdata.all)[mdata.all==1];
  c(substring(sel.nms, 0, nchar(sel.nms)-4), colnames(meta.mat));
}

# single.type return logFC or p value for individual data analysis
GetMetaResultMatrix<-function(single.type="fc"){
  if(single.type == "fc"){
    meta.mat <- cbind(fc.mat, meta.mat);
  }else{
    meta.mat <- cbind(pval.mat, meta.mat);
  }
  # display at most 500 genes
  if(nrow(meta.mat) > 500){
    meta.mat <- meta.mat[1:500,]; # already sorted based on meta-p values
  }
  meta.mat <-signif(as.matrix(meta.mat), 5);
  meta.mat;
}

GetMetaStat<-function(){
  return (meta.stat$stat);
}

GetMetaStatNames<-function(){
  return (names(meta.stat$stat));
}

combinePvals <- function(pvalonesided,nrep,BHth=0.05, method) {
  listres=vector("list",3);
  nbstudies=length(pvalonesided);
  nbreptot=sum(nrep);
  if (nbreptot <2) {
    stop("Error: the argument \"nrep\" must be a vector with at least two values higher than 1")
  } 
  
  weight=sqrt(nrep/nbreptot);
  fstatmeta=function(g){
    vecptime=unlist(lapply(pvalonesided, FUN = function(x) x[g]));
    vec = qnorm(1 - vecptime);
    stattestg = sum(weight[1:length(pvalonesided)] * vec[1:length(pvalonesided)], na.rm = TRUE);
    stattestg;
  }
  
  fishersum <- function(pvec){
    return(sum(-2*log(pvec)))
  }
  
  if(method=="stouffer"){
    statpvalc=-unlist(lapply(rep(1:length(as.vector(pvalonesided[[1]])), 1), function(x) fstatmeta(x)));
    rpvalpvalc=2*(1-pnorm(abs(statpvalc)));
  }else{ # fisher
    data <- data.frame(pvalonesided);
    #data[data == 0] <- 1*10^-10;
    
    #note, p value are calculated for one side
    # pt (lower.tail=T by default) which tests if group A < group B
    # for one side
    fsum1 <- apply(data, 1, fishersum);
    rpvalpvalc1 = 1-pchisq(fsum1, df=(ncol(data)*2));
    
    # for the other side
    data <- 1-data;
    fsum2 <- apply(data, 1, fishersum);
    rpvalpvalc2 = 1-pchisq(fsum2, df=(ncol(data)*2));
    
    # report the min of two direction calculation
    rpvalpvalc <- pmin(rpvalpvalc1, rpvalpvalc2);
    
    # adding direction information
    statpvalc <- pmax(fsum1, fsum2);
    # if A<B sig, then it should be negative 
    statpvalc[statpvalc == fsum1]<- -statpvalc[statpvalc == fsum1];
  }
  
  rpvalpvalc <- p.adjust(rpvalpvalc,method="BH");
  res=which(rpvalpvalc<=BHth);
  listres[[1]]=res
  listres[[2]]=statpvalc;
  listres[[3]]=rpvalpvalc
  names(listres)=c("DEindices", "CombinedStat", "CombinedP")
  listres
}

#combine effect size
combineES <- function (ES, varES, BHth = 0.05, method){
  if(method == "rem"){
    useREM = TRUE;
  }else{
    useREM = FALSE;
  }
  
  num.studies <- dim(ES)[2];
  
  Qvals <- f.Q.NA(ES, varES)
  if (useREM) {
    varES <- varES + tau2.NA(Qvals, num.studies, my.weights = 1/varES)
  }
  wt <- 1/varES
  MUvals <- rowSums(ES * wt, na.rm = TRUE)/rowSums(wt, na.rm = TRUE)
  MUsES <- sqrt(abs(1/rowSums(wt, na.rm = TRUE)))
  zSco <- MUvals/MUsES
  rpvalESc = 2 * (1 - pnorm(abs(zSco)))
  res = which(p.adjust(rpvalESc, method = "BH") <= BHth);
  listres <- list();
  listres[[1]] = res
  listres[[2]] = zSco
  listres[[3]] = MUvals; # pool effect size
  listres[[4]] = wt; # wt for each studies, it is matrix with one column for each studies
  names(listres) = c("DEindices", "TestStatistic", "PooledEffectSize", "Weights")
  listres
}

# first two biased ES, var, last two unbiased

f.Q.NA = function(dadj, varadj) {
  w <- 1/varadj
  tmp1 <- w * dadj
  mu <- rowSums(tmp1, na.rm = TRUE)/rowSums(w, na.rm = TRUE)
  Q <- rowSums(w * (dadj - mu)^2, na.rm = TRUE)
}

tau2.NA <- function(Q, num.studies, my.weights) {
  vwts <- rowSums(my.weights, na.rm = TRUE)
  tmp2 <- rowSums(my.weights^2, na.rm = TRUE)
  tau2 <- pmax(0, (Q - (num.studies - 1))/(vwts - tmp2/vwts))
  return(tau2)
}

# prepare data for heatmap plotting include
# 1. meta info (top bars)
# 2. expression matrix
# 3. function annotation (left bars)
# all matrix will be combined, 
# 1 and 2 separated by a row of 'null' 
# 3 and 1+2 separated by a column of 'null'
PrepareMetaHeatmapJSON <- function(){
  gene.vec <- rownames(meta.mat);
  datanm.vec <- names(mdata.all)[mdata.all==1];
  
  inmex.meta <- readRDS("inmex_meta.rds");
  dat.inx <- inmex.meta$data.lbl %in% datanm.vec;
  dat <- inmex.meta$plot.data[gene.vec, dat.inx, drop=F]; 
  
  # scale each gene for each dataset
  dat <- t(scale(t(dat)));
  
  # now need to remove na or constant rows
  dat <- na.omit(dat);
  # check for columns with all constant (var =0)
  varCol <- apply(dat, 1, var, na.rm=T);
  constCol <- (varCol == 0 | is.na(varCol));
  dat <- dat[!constCol, ];
  
  anot.res <- list();
  ids <- rownames(dat);
  if(inmex.meta$id.type == "entrez"){   
    anot.res <- doEntrezIDAnot(ids);
  }else{ # no annotation, then use the default feature ID
    anot.res$gene_id <- anot.res$symbol <- anot.res$name <- ids; 
  }
  
  data.nms <- as.factor(inmex.meta$data.lbl[dat.inx]);
  cls.lbls <- as.factor(inmex.meta$cls.lbl[dat.inx]);
  
  # setup annotation info
  data.nms <- as.character(data.nms);
  datasets <- substr(as.character(data.nms), 0, nchar(data.nms)-4);
  
  sel.nms <- names(mdata.all)[mdata.all==1];
  if(length(sel.nms) > 1){
    annotation <- data.frame(class= cls.lbls, dataset = as.factor(datasets));
  }else{ # single data
    annotation <- data.frame(class= cls.lbls);
  }
  rownames(annotation) <- colnames(dat);
  
  ####
  sig.ids <- rownames(dat);
  if(inmex.method != "votecount"){
    stat.pvals <- as.numeric(meta.mat[,2]);
    stat.fc = as.numeric(meta.mat[,1]);
  }else{
    stat.pvals <- as.numeric(meta.mat[,1]);
    stat.fc = as.numeric(meta.mat[,1]);
  }
  
  orig.smpl.nms <- colnames(dat);
  orig.gene.nms <- rownames(dat);
  
  # do clustering and save cluster info
  if(nrow(dat)> 1){
    dat.dist <- dist(dat);
    gene.ward.ord <- hclust(dat.dist, "ward.D")$order;
    gene.ward.rk <- match(orig.gene.nms, orig.gene.nms[gene.ward.ord]);
    gene.ave.ord <- hclust(dat.dist, "ave")$order;
    gene.ave.rk <- match(orig.gene.nms, orig.gene.nms[gene.ave.ord]);
    gene.single.ord <- hclust(dat.dist, "single")$order;
    gene.single.rk <- match(orig.gene.nms, orig.gene.nms[gene.single.ord]);
    gene.complete.ord <- hclust(dat.dist, "complete")$order;
    gene.complete.rk <- match(orig.gene.nms, orig.gene.nms[gene.complete.ord]);
    
    dat.dist <- dist(t(dat));
    smpl.ward.ord <- hclust(dat.dist, "ward.D")$order;
    smpl.ward.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ward.ord])
    smpl.ave.ord <- hclust(dat.dist, "ave")$order;
    smpl.ave.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ave.ord])
    smpl.single.ord <- hclust(dat.dist, "single")$order;
    smpl.single.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.single.ord])
    smpl.complete.ord <- hclust(dat.dist, "complete")$order;
    smpl.complete.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.complete.ord])
  }else{
    # force not to be single element vector which will be scaler
    stat.pvals <- matrix(stat.pvals);
    gene.ward.rk <- gene.ave.rk <- gene.single.rk <- gene.complete.rk <- matrix(1);
    smpl.ward.rk <- smpl.ave.rk <- smpl.single.rk <- smpl.complete.rk <- 1:ncol(dat);
  }
  
  gene.cluster <- list(
    pval = stat.pvals,
    fc = stat.fc, 
    ward = gene.ward.rk,
    average = gene.ave.rk,
    single = gene.single.rk,
    complete = gene.complete.rk
  );
  
  sample.cluster <- list(
    ward = smpl.ward.rk,
    average = smpl.ave.rk,
    single = smpl.single.rk,
    complete = smpl.complete.rk
  );
  
  # prepare meta info    
  # 1) convert meta.data info numbers
  # 2) match number to string (factor level)
  meta <- annotation;
  grps <- colnames(meta);
  nmeta <- meta.vec <- NULL;
  uniq.num <- 0;
  for (i in 1:ncol(meta)){
    cls <- meta[,i];
    grp.nm <- grps[i];
    meta.vec <- c(meta.vec, as.character(cls))
    # make sure each label are unqiue across multiple meta data
    ncls <- paste(grp.nm, as.numeric(cls)); # note, here to retain ordered factor
    nmeta <- c(nmeta, ncls);
  }
  
  # convert back to numeric 
  nmeta <- as.numeric(as.factor(nmeta))+99;
  unik.inx <- !duplicated(nmeta)   
  
  # get corresponding names
  meta_anot <- meta.vec[unik.inx]; 
  names(meta_anot) <- nmeta[unik.inx]; # name annotatation by their numbers
  
  nmeta <- matrix(nmeta, ncol=ncol(meta), byrow=F);
  colnames(nmeta) <- grps;
  
  # for each gene/row, first normalize and then tranform real values to 30 breaks 
  res <- t(apply(dat, 1, function(x){as.numeric(cut(x, breaks=30))}));
  
  # note, use {} will lose order; use [[],[]] to retain the order
  # single element vector will be converted to scalar, not array, need to prevent that
  gene.id = anot.res$symbol; if(length(gene.id) ==1) { gene.id <- matrix(gene.id) };
  gene.entrez = anot.res$gene_id; if(length(gene.entrez) ==1) { gene.entrez <- matrix(gene.entrez) };        
  gene.name = anot.res$name; if(length(gene.name) ==1) { gene.name <- matrix(gene.name) };
  
  json.res <- list(
    data.type = "array",
    gene.id = as.character(anot.res$symbol),
    gene.entrez = gene.entrez,
    gene.name = anot.res$name,
    gene.cluster = gene.cluster,
    sample.cluster = sample.cluster,
    sample.names = orig.smpl.nms,
    meta = data.frame(nmeta),
    meta.anot = meta_anot,
    data.lbl = inmex.meta$data.lbl,
    data = res
  );
  
  return(json.res);
}


##################################
# functions for estimating Cochran’s Q
##################################


#computes Cochran’s Q gene by gene
#dadj and varadj must be matrices, in which every study is a column,
#every row a gene
f.Q <- function(dadj, varadj){
  w<-1/varadj
  tmp1<-w*dadj
  mu<-rowSums(tmp1)/rowSums(w)
  Q<-rowSums(w*(dadj - mu)^2)
}

qc.metaDensity<- function(imgNm, dpi=72, format, factor){
  library("ggplot2")
  inmex.meta <- readRDS("inmex_meta.rds");
  dat = inmex.meta$data
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  dpi = as.numeric(dpi)
  
  df = data.frame(inmex.meta$data, stringsAsFactors = FALSE)
  df = stack(df)
  if(factor != "NA"){
    #Factor = as.vector(dataSet$meta.info[,factor])
  }else{
    #Factor = dataSet$meta.info[,1];
  }
  Factor=inmex.meta$data.lbl
  
  conv = data.frame(ind=colnames(inmex.meta$data), class=Factor)
  conv$ind=gsub("-", ".", conv$ind)
  df1 = merge(df, conv, by="ind")
  Cairo(file=imgNm, width=10, height=6, type=format, bg="white", dpi=dpi, unit="in");
  g =ggplot(df1, aes(x=values)) + geom_line(aes(color=class, group=ind), stat="density", alpha=0.3) + geom_line(aes(color=class), stat="density", alpha=0.6, size=1.5)
  print(g)
  dev.off();
}

effectsize <- function(tstat,ntilde,m){
  cm=gamma(m/2)/(sqrt(m/2)*gamma((m-1)/2))
  d=tstat/sqrt(ntilde)
  dprime=cm*d
  terme1=m/((m-2)*ntilde)
  vard=terme1+d^2*(terme1*ntilde-1/cm^2)
  vardprime=cm^2*(terme1+dprime^2*(terme1*ntilde-1/cm^2))
  result=cbind(d,vard,dprime,vardprime)
  colnames(result)=c("d","vard","dprime","vardprime")
  result
}

#For internal GSEA of meta-analysis of gene sets

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

PlotMetaPhm <-function(cmpdNm, dpi=72){
  allmat = readRDS("allmat.rds");
  
  fileNm = paste("Path_", cmpdNm, ".png", sep="");
  current.geneset <- readRDS("current_geneset.rds")
  ids=current.geneset[[cmpdNm]];
  subset = allmat[which(rownames(allmat) %in% ids),]
  if(length(subset)<1){
    subset = allmat[which(rownames(allmat) %in% ids),]
  }
  library(RColorBrewer);
  library(pheatmap)
  dims = dim(subset)
  rnms = rownames(subset)
  cnms = colnames(subset) 
  m <- mapply(subset, FUN=as.numeric)
  subset <- matrix(data=m, ncol=dims[2], nrow=dims[1])
  rownames(subset) = rnms
  colnames(subset) = cnms
  if(inmex.method %in% c("effectsize", "merge")){
    subset = subset[complete.cases(subset), ];
  }
  
  subset[is.na(subset)] <- 0
  my_palette <- colorRampPalette(c("green", "black", "red"))(n = 30)
  my_palette = c("#d3d3d3", my_palette)
  bk2 = unique(c(seq(0,0.9999, length=2), 1, seq(2,30, length=29)));
  inmex.meta <- readRDS("inmex_meta.rds");
  ann = data.frame(Class=inmex.meta$cls.lbl, Dataset=inmex.meta$data.lbl)
  rownames(ann) = colnames(subset)
  Cairo(file=fileNm, width=800, height=700, type="png", bg="white",unit="px",dpi=72)
  hm = pheatmap(subset, color = my_palette, breaks = bk2, annotation_col=ann, show_rownames = FALSE, main = cmpdNm, border_color=NA)
  print(hm)
  dev.off()
}

# for gene set-level meta-analysis

SetupMetaGSEAStats <- function(name, netNm, BHth, mType, curr.geneset, lib){
  
  inmex.de <- list();
  allmat <- readRDS("allMeta.mat.rds");
  allmat.vec <- rownames(allmat);
  meta.mat = metaset.mat
  metade.genes <- rownames(meta.mat);
  pval.mat <- fc.mat <- matrix(nrow=nrow(meta.mat), ncol=sum(mdata.all==1));
  fc.list <- split(rep(" ", length(allmat.vec)), allmat.vec);
  
  
  current.geneset = curr.geneset[!duplicated(names(curr.geneset))]
  inx =names(current.geneset) %in% rownames(meta.mat)  ;
  
  resTable = meta.mat
  current.mset = current.geneset[inx];
  
  inmex.meta <- readRDS("inmex_meta.rds");
  
  ora.vec <- rownames(inmex.meta$data)
  ora.nms <- doEntrez2SymbolMapping(rownames(inmex.meta$data))
  
  hits.query <- lapply(current.mset, 
                       function(x) {
                         unique(ora.nms[ora.vec%in%unlist(x)]);
                       });  
  saveRDS(hits.query, "hits_query.rds");
  
  set.num = unlist(lapply(current.mset, function(x){length(unique(x))}), use.names=TRUE);
  names(hits.query) <- names(current.mset);
  hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=TRUE);
  
  vote.bool = "false"
  meta.matcolinx = 2;
  enr.score = "NA"   
  
  padj <- p.adjust(as.vector(meta.mat[,meta.matcolinx]),method="BH");
  if(mType == "network"){
    json.res <- list(
      fun.anot = hits.query,
      fun.ids = as.vector(rownames(meta.mat)),
      fun.pval = as.vector(meta.mat[,meta.matcolinx]),
      fun.padj = padj,
      hit.num = hit.num,
      total= set.num
    );
  }else{
    json.res <- list(
      hits = hit.num,
      total= set.num,
      enr.pval= as.vector(meta.mat[,meta.matcolinx]),
      enr.padj= padj,
      enr.names= as.vector(rownames(meta.mat)),
      cls.lbl=inmex.meta$cls.lbl,
      smps.lbl=smps.vec,
      data.lbl = inmex.meta$data.lbl,
      path.lbl = rownames(meta.mat),
      enr.score = as.vector(meta.mat[,1]),
      isVote = vote.bool
    );
  }
  
  
  library(RJSONIO);
  json.mat <- RJSONIO::toJSON(json.res, .na='null', pretty=TRUE);
  json.nm <- paste0(name, ".json");
  
  sink(json.nm);
  cat(json.mat);
  sink();
  if(mType == "network"){
    res.mat<-matrix(0, nrow=length(names(current.mset)), ncol=4);
    rownames(res.mat)<-names(current.mset);
    colnames(res.mat)<-c("Total","Hits", "P.Value", "FDR");
    res.mat[,"Total"] = set.num
    res.mat[,"Hits"] = hit.num;
    res.mat[,"P.Value"] = meta.mat[,meta.matcolinx];
    res.mat[,"FDR"] = padj;
    enr.mat <<- res.mat
    res <- data.frame(Name=as.vector(rownames(meta.mat)), Total=set.num, Hits= hit.num, EnrichmentScore=as.vector(meta.mat[,1]), Pval=as.vector(meta.mat[,meta.matcolinx]), Padj = padj);
    list.genes <<- allmat.vec
    SetListNms();
    netnm <- paste0(netNm, ".json");
    PrepareEnrichNet(netNm, "meta", "mixed");
  }else{
    res.mat<-matrix(0, nrow=length(names(current.mset)), ncol=5);
    colnames(res.mat)<-c("Name", "Total","Hits", "P.Value", "FDR");
    res.mat[,"Name"] = names(current.mset);
    res.mat[,"Total"] = set.num
    res.mat[,"Hits"] = hit.num;
    res.mat[,"P.Value"] = meta.mat[,meta.matcolinx];
    res.mat[,"FDR"] = padj;
    write.csv(res.mat, file=paste("meta_sig_genesets_", lib, ".csv", sep=""), row.names=F);
  }
  return(1)
}

CalculateGsNet <- function(name, netNm, type, mType, db){
  res = PerformMetaPathCombine(name, netNm, "pval", db, mType,0.05)
  return(1)
}
