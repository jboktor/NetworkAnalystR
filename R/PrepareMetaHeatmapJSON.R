# prepare data for heatmap plotting include
# 1. meta info (top bars)
# 2. expression matrix
# 3. function annotation (left bars)
# all matrix will be combined, 
# 1 and 2 separated by a row of 'null' 
# 3 and 1+2 separated by a column of 'null'
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PrepareMetaHeatmapJSON
#' @export 
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
