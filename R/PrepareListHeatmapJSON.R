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
#' @rdname PrepareListHeatmapJSON
#' @export 
PrepareListHeatmapJSON <- function(){
  sig.ids <- rownames(dataSet$prot.mat);
  gene.symbols=doEntrez2SymbolMapping(sig.ids)
  stat.pvals <- dataSet$prot.mat[,1]
  
  expval <- 0
  expval <- sum(dataSet$prot.mat)
  
  # scale each gene 
  #data.stat <- readRDS("data.stat");
  #hit.inz = sig.ids %in% rownames(data.stat);
  #sig.ids = sig.ids[hit.inz];
  dat <- dataSet$prot.mat
  
  # now pearson and euclidean will be the same after scaleing
  dat.dist <- dist(dat); 
  
  orig.smpl.nms <- colnames(dat);
  orig.gene.nms <- rownames(dat);
  
  # prepare meta info    
  # 1) convert meta.data info numbers
  # 2) match number to string (factor level)
  
  grps <- "datalist1"
  cls <- "datalist1"
  
  # convert back to numeric 
  
  # for each gene/row, first normalize and then tranform real values to 30 breaks
  if(expval !=0){
    dat_pos = as.matrix(dat[sign(dat[,1]) == 1,])
    dat_neg = as.matrix(dat[sign(dat[,1]) == -1,])
    if(nrow(dat_pos) == 0){
      res <- apply(unname(dat), 2, function(x){
        y =log(abs(x)) + 0.000001
        16-as.numeric(cut(y, breaks=15))
      });
    }else if(nrow(dat_neg) == 0){
      res <- apply(unname(dat), 2, function(x){
        y =log(x) + 0.000001
        15+as.numeric(cut(y, breaks=15))
      });
    }else{
      res_pos <- apply(unname(dat_pos), 2, function(x){
        y =log(x) + 0.000001
        as.numeric(cut(y, breaks=15))+15
      });
      res_neg <- apply(unname(dat_neg), 2, function(x){
        y =log(abs(x)) + 0.000001
        16 - as.numeric(cut(y, breaks=15))
      });
      res = rbind(res_pos, res_neg);
    }
  }else{
    zero.inx = dataSet$prot.mat == 0
    res = dataSet$prot.mat;
    res[zero.inx] = 32
  }
  
  res_list <- list()
  for(i in 1:length(res)){
    res_list[[i]] <- list(res[i])
  }
  
  # note, use {} will lose order; use [[],[]] to retain the order
  
  nmeta = list(100)
  nmeta.anot = list()
  
  nmeta.anot["datalist1"] = nmeta[1]
  
  nmeta = list(nmeta)
  names(nmeta) = "datalists"
  
  json.res <- list(
    data.type = "singlelist", 
    gene.id = gene.symbols,
    gene.entrez = sig.ids,
    gene.name = gene.symbols,
    gene.cluster = 1,
    sample.cluster = 1,
    sample.names = list("datalist1"),
    meta = nmeta,
    meta.anot = nmeta.anot,
    data = res_list,
    expval = expval
  );
  rownames(dat) = gene.symbols
  write.csv(dat,"heatmap_matrix.csv", row.names=TRUE)
  return(json.res);
}
