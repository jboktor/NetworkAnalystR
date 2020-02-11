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
#' @rdname PrepareMultiListHeatmapJSON
#' @export 
PrepareMultiListHeatmapJSON <- function(){
  sel.nms <- names(mdata.all)
  expval<-0;
  for(i in 1:length(sel.nms)){
    dataNm <- sel.nms[i];
    dataSet <- readRDS(dataNm);
    len <- nrow(dataSet$prot.mat)
    if(i == 1){
      expval <- sum(dataSet$prot.mat)
      gene_list <-rownames(dataSet$prot.mat)
    }else{
      gene_list <-c(gene_list, rownames(dataSet$prot.mat))
      expval <- expval + sum(dataSet$prot.mat)
    }
  }
  
  gene_list <- unique(gene_list)
  allmat = matrix(NA, nrow=length(gene_list), ncol=length(sel.nms))
  rownames(allmat) = gene_list
  
  for(i in 1:length(sel.nms)){
    dataName <- sel.nms[i];
    dataSet <- readRDS(dataName);
    cols <- colnames(allmat)[colnames(allmat) %in% dataName]
    if(expval ==0){
      rows <- which(rownames(allmat) %in% rownames(dataSet$prot.mat))
      inx <-match(rownames(allmat) ,rownames(dataSet$prot.mat))
      allmat[, i] <- as.vector(dataSet$prot.mat)[inx]
    }else{
      rows <- which(rownames(allmat) %in% rownames(dataSet$prot.mat))
      inx <-match(rownames(allmat) ,rownames(dataSet$prot.mat))
      allmat[, i] <- as.vector(dataSet$prot.mat)[inx]
    }
  } 
  colnames(allmat) = sel.nms 
  inx <- apply(allmat, 1, function(x){sum(is.na(x))});  
  ord.inx <- order(inx)
  allmat = allmat[ord.inx,]
  gene.symbols = doEntrez2SymbolMapping(rownames(allmat))
  
  na.inx = is.na(allmat)
  zero.inx = allmat == 0
  
  allmatb = allmat
  
  allmatb[na.inx]=0
  allmatb[zero.inx]=1
  rownames(allmatb) = gene.symbols
  write.csv(allmatb,"heatmap.csv", row.names=TRUE)  
  
  if(expval != 0){
    pos.inx = allmat>0 & !na.inx
    neg.inx = allmat<0 & !na.inx
    allmat[neg.inx] = 16 - as.numeric(cut(log(abs(allmat[neg.inx])) , breaks=15))
    allmat[pos.inx] = 15 + as.numeric(cut(log(allmat[pos.inx]) , breaks=15))
    allmat[zero.inx] = 32
  }else{
    zer.inx = allmat == 0 & !na.inx
    nb <- apply(allmat, 1, function(x){sum(!is.na(x))});
    for(i in 1:nrow(allmat)){
      row = allmat[i,]
      inx = row == 0
      allmat[i, inx] = nb[i]
    }
    allmat[zer.inx] <- try(15 + as.numeric(cut(allmat[zer.inx] , breaks=15)));
    if(class(allmat[zer.inx]) == "try-error") {
      allmat[zer.inx] = 32
    }else{
      15 + as.numeric(cut(allmat[zer.inx] , breaks=15))
    }
  }
  allmat[na.inx] = 31;
  res_list <- list();
  for(i in 1:nrow(allmat)){
    res_list[i] <- list(unname(allmat[i,]))
  }
  
  nmeta = as.numeric(as.factor(colnames(allmat))) + 99
  nmeta.anot = list()
  
  for(i in 1:length(unique(nmeta))){
    nmeta.anot[[colnames(allmat)[i]]] = nmeta[i]
  }
  nmeta = list(nmeta)
  names(nmeta) = "datalists"
  
  json.res <- list(
    data.type = "mutlilist",
    gene.id = gene.symbols,
    gene.entrez = rownames(allmat),
    gene.name = rownames(allmat),
    gene.cluster = 1,
    sample.cluster = 1,
    sample.names = colnames(allmat),
    meta = nmeta,
    meta.anot = nmeta.anot,
    data.lbl = "NA",
    data = res_list,
    expval = expval
  );
  
  return(json.res);
}
