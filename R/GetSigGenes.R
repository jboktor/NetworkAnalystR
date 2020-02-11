# update result based on new cutoff
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param res.nm PARAM_DESCRIPTION
#' @param p.lvl PARAM_DESCRIPTION
#' @param fc.lvl PARAM_DESCRIPTION
#' @param update PARAM_DESCRIPTION, Default: T
#' @param inx PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetSigGenes
#' @export 
GetSigGenes <-function(res.nm, p.lvl, fc.lvl, update=T, inx){
  total = nrow(dataSet$resTable);
  resTable <- dataSet$resTable;
  filename <- dataSet$filename;
  filename <- paste(filename, "_", res.nm, ".csv", sep="");
  if (update){
    current.msg <<- "";
  }
  # select based on p-value
  if(dataSet$type == "array"){
    hit.inx.p <- resTable$adj.P.Val <= p.lvl; 
  } else {
    hit.inx.p <- resTable$adj.P.Val <= p.lvl; 
  }
  
  resTable<-resTable[hit.inx.p,,drop=F];
  if (nrow(resTable) == 0){
    current.msg <<- paste(current.msg, "No significant genes were identified using the given design and cutoff."); 
  }
  # now rank by logFC, note, the logFC for each comparisons 
  # are returned in resTable before the AveExpr columns 
  # for two-class, only one column, multiple columns can be involved
  # for > comparisons - in this case, use the largest logFC among all comparisons
  #if (fc.lvl > 0){ # further filter by logFC
  if (dataSet$de.method=="limma"){
    hit.inx <- which(colnames(resTable) == "AveExpr");
  } else if (dataSet$de.method=="deseq2"){
    hit.inx <- which(colnames(resTable) == "baseMean");  
  } else {
    hit.inx <- which(colnames(resTable) == "logCPM");
  }
  maxFC.inx <- hit.inx - 1; # not sure if this is also true for edgeR
  logfc.mat <- resTable[,1:maxFC.inx, drop=F];
  pos.mat <- abs(logfc.mat);
  fc.vec <- apply(pos.mat, 1, max);
  hit.inx.fc <- fc.vec >= fc.lvl;
  resTable<-resTable[hit.inx.fc,,drop=F];
  if (nrow(resTable) == 0){
    current.msg <<- paste(current.msg, "No significant genes were identified using the given design and cutoff."); 
    
  }
  #}
  
  ### Note, rowname of resTable must be entrez ID
  
  de.Num <- nrow(resTable);
  
  # display at most 5000 genes for the server (two main reasons)
  # 1) should not have more 22% (human: 23000) DE of all genes (biological)
  # 2) IE canvas can display no more than 6800 pixels (computational)
  if (nrow(resTable) > 5000){
    resTable <- resTable[1:5000,];
    current.msg <<- paste(current.msg, " Due to computational constraints, only top 5000 genes will be used. ", collapse="\n");
  }
  
  # may need to update data, class and meta.info
  data <- dataSet$data.norm;
  cls <- dataSet$cls; 
  meta.info <- dataSet$meta.info;
  grp.nms <- levels(cls);
  
  hit.inx <- cls %in% grp.nms;
  if (sum(hit.inx) < length(hit.inx)){
    current.msg <<- paste(current.msg, "Only groups selected for comparisons: ", paste(grp.nms, collapse=", "), "are included.");
    cls <- factor(cls[hit.inx]);
    cls.lvls <- levels(cls);
    data <- data[,hit.inx];
    meta.info <- dataSet$meta.info[hit.inx,];
  }
  
  saveRDS(data, file="data.stat");
  dataSet$resTable = dataSet$resTable[order(dataSet$resTable$adj.P.Val),] 
  dataSet$resTable = dataSet$resTable[which(!rownames(dataSet$resTable) %in% rownames(resTable)),]
  dataSet$resTable = rbind(resTable, dataSet$resTable);
  
  dataSet$sig.mat <- resTable;
  if (dataSet$annotated){ # annotated to entrez
    anot.id <- rownames(dataSet$resTable);
    gene.anot <- doEntrezIDAnot(anot.id);
    write.csv(cbind(EntrezID=anot.id, signif (dataSet$resTable,5), Symbols = gene.anot$symbol,  Name=gene.anot$name), row.names=F, file=filename);
  } else if (file.exists("annotation.rds")){ # annotation information available
    anot.id <- readRDS("annotation.rds");
    feature.vec <- rownames(dataSet$resTable);
    entrez.vec <- anot.id[feature.vec];
    gene.anot <- doEntrezIDAnot(entrez.vec);
    write.csv(cbind(signif (dataSet$resTable,5), EntrezID=entrez.vec, Symbols = gene.anot$symbol,  Name=gene.anot$name), row.names=F, file=filename);
    rownames(gene.anot) <- feature.vec;
  } else {
    gene.anot <- NULL;
    write.csv(signif(resTable,5), file=filename);
  }
  if(is.null(gene.anot)){
    dataSet$sig.genes.symbols <- rep("NA",nrow(resTable));
  }else{
    dataSet$sig.genes.symbols <- gene.anot$symbol;
  }
  dataSet$cls.stat <- cls;
  dataSet$meta.stat <- meta.info;
  
  # now do protein mapping for network only applicable for annotated
  
  dataSet$name <- res.nm;
  
  gene <- rownames(resTable);
  
  logFC <- unname(logfc.mat[,1]);
  geneList <- paste(gene, logFC, collapse="\n");
  up = nrow(resTable[which(logfc.mat[,selectedFactorInx]> fc.lvl),])
  down = nrow(resTable[which(logfc.mat[,selectedFactorInx]< -fc.lvl),])
  
  dataSet <<- dataSet;
  data.norm <- dataSet$data.norm
  colnames(data.norm) = NULL
  lst = list(colnames(dataSet$data.norm),data.norm, dataSet$meta.info, dataSet$resTable, rownames(data.norm), org=data.org)
  require(RJSONIO)
  json.obj <- toJSON(lst);
  sink("NetworkAnalyst_matrix.json");
  cat(json.obj);
  return(c(filename, de.Num, geneList, total, up, down));
}
