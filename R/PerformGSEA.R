#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param file.nm PARAM_DESCRIPTION
#' @param fun.type PARAM_DESCRIPTION
#' @param netNm PARAM_DESCRIPTION
#' @param mType PARAM_DESCRIPTION
#' @param selectedFactorInx PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformGSEA
#' @export 
PerformGSEA<- function(file.nm, fun.type, netNm, mType, selectedFactorInx){
  
  
  if(tolower(fun.type) == 'kegg'){ 
    current.geneset <- LoadKEGGLib();
  }else if(tolower(fun.type) == 'reactome'){ 
    current.geneset <- LoadREACTOMELib();
  }else if(tolower(fun.type) == 'motif'){ 
    current.geneset <- LoadMotifLib();
  }else{ # GO
    current.geneset <- LoadGOLib(fun.type);
  }
  
  require("fgsea");
  
  if(anal.type == "onedata"){
    datnorm = dataSet$data.norm
    sampleNms = colnames(dataSet$data.norm);
    rankedVec<- ComputeRankedVec(dataSet, selectedFactorInx);
    
  }else{
    if(dataSet$name != selDataNm){
      dataSet <- readRDS(selDataNm);
    }
    datnorm=dataSet$data
    sampleNms = colnames(dataSet$data);
    ds = inmex.ind[selDataNm][[1]]
    rankedVec <- ComputeRankedVec(dataSet, 1);
  }
  
  fgseaRes <- fgsea(pathways = current.geneset, 
                    stats = rankedVec,
                    minSize=15,
                    maxSize=500,
                    nperm=10000);
  
  fgseaRes <- fgseaRes[!duplicated(fgseaRes$pathway),]
  
  rownames(fgseaRes) = make.names(fgseaRes$pathway, unique=TRUE)
  fgseaRes = fgseaRes[,c("size","ES", "pval", "pathway", "padj")]
  
  if(nrow(fgseaRes)<1){
    require(RJSONIO);
    SetListNms()
    initsbls <- doEntrez2SymbolMapping(list.genes);
    names(initsbls) <- list.genes
    netData <- list(sizes=listSizes, genelist=initsbls);
    netName <- paste0(netNm, ".json");
    sink(netName);
    cat(toJSON(netData));
    sink();
    return(0);
  }
  
  fgseaRes <- fgseaRes[order(fgseaRes$pval),]
  fgseaRes <<- fgseaRes
  if(nrow(fgseaRes[which(fgseaRes$pval < 0.05),])<20 ){
    if(nrow(fgseaRes)>20){
      fgseaRes = fgseaRes[c(1:20),]
    }
  }else{
    fgseaRes[which(fgseaRes$pval < 0.05),]
  } 
  
  inx <- which(names(current.geneset) %in% fgseaRes$pathway);
  current.mset = current.geneset[inx]
  current.mset = current.mset[!duplicated(names(current.mset))]
  
  ora.vec <- names(rankedVec)
  ora.nms <- doEntrez2SymbolMapping(ora.vec)
  
  hits.query <- lapply(current.mset, 
                       function(x) {
                         unique(ora.nms[ora.vec%in%unlist(x)]);
                       });
  saveRDS(hits.query, "hits_query.rds");
  set.num = unlist(lapply(current.mset, function(x){length(unique(x))}), use.names=TRUE);
  names(hits.query) <- names(current.mset);
  hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=TRUE);
  fgseaRes$hits = hit.num[which(fgseaRes$pathway  %in% names(hit.num))] 
  fgseaRes$total = set.num[which(fgseaRes$pathway %in% names(set.num))]
  
  fgseaRes = fgseaRes[which(fgseaRes$hits>1),]
  fgseaRes = fgseaRes[which(fgseaRes$hits<500),]
  fgseaRes = fgseaRes[which(fgseaRes$total<2000),]
  if(nrow(fgseaRes)<1){
    require(RJSONIO);
    SetListNms();
    initsbls = doEntrez2SymbolMapping(list.genes)
    names(initsbls) = list.genes
    netData <- list(sizes=listSizes, genelist=initsbls);
    netName = paste0(netNm, ".json");
    sink(netName);
    cat(toJSON(netData));
    sink();
    return(0);
  }
  
  fgseaRes=fgseaRes[order(fgseaRes$pval),]
  if(nrow(fgseaRes[which(fgseaRes$pval < 0.05),])<20 ){
    if(nrow(fgseaRes)>20){
      fgseaRes = fgseaRes[c(1:20),]
    }
  }else{
    fgseaRes = fgseaRes[which(fgseaRes$padj < 0.05),]
  } 
  
  fgseaRes <- data.frame(fgseaRes, stringsAsFactors=FALSE)
  
  #get gene symbols
  current.msg <<- "Functional enrichment analysis was completed";
  
  # write json
  require(RJSONIO);
  fun.anot = hits.query; 
  fun.pval = fgseaRes[,3]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  #fun.pval<-signif(fun.pval,5);  
  fun.padj = fgseaRes[,5]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  #fun.padj<-signif(fun.padj,5);  
  es.num = fgseaRes[,2]; if(length(es.num) ==1) { es.num <- matrix(es.num) };
  fun.ids <- as.vector(current.setids[names(fun.anot)]); 
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  
  json.res <- list(
    fun.link = current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    pathname = fgseaRes[,"pathway"],
    es.num = es.num,
    hits = fgseaRes[,"hits"],
    total = fgseaRes[,"total"],
    cls = dataSet$meta.info,
    sample.nms = sampleNms       
  );
  
  if(mType == "network"){
    res.mat<-matrix(0, nrow=length(fun.pval), ncol=4);
    colnames(res.mat)<-c("Total","Hits", "P.Value", "FDR");
    res.mat[,"Total"] = fgseaRes[,"total"];
    res.mat[,"Hits"] = fgseaRes[,"hits"];
    res.mat[,"P.Value"] = fgseaRes[,"pval"];
    res.mat[,"FDR"] = fgseaRes[,"padj"];
    res.mat = data.matrix(data.frame(res.mat, stringsAsFactors=FALSE));
    rownames(res.mat) = fgseaRes[,"pathway"];
    enr.mat <<- res.mat;
    list.genes <<- doEntrez2SymbolMapping(rownames(dataSet$sig.mat));
    SetListNms();
    PrepareEnrichNet(netNm, "meta", "mixed");
    file.nm <- gsub("gsea", "enrichment", file.nm)
    json.mat <- toJSON(json.res, .na='null');
    json.nm <- paste(file.nm, ".json", sep="");
  }else{
    
    json.mat <- toJSON(json.res, .na='null');
    json.nm <- paste(file.nm, ".json", sep="");
  }
  
  sink(json.nm)
  cat(json.mat);
  sink();
  
  # write csv
  
  fgseaRes <<- fgseaRes 
  ftype = fun.type
  if(fun.type %in% c("bp", "mf", "cc")){
    ftype = paste0("go_", fun.type);
  }
  csvDf <- data.frame(Name=fgseaRes$pathway, Total=fgseaRes$total, Hits=fgseaRes$hits, EnrichmentScore=fgseaRes$ES, Pval=fgseaRes$pval, Padj=fgseaRes$padj);
  write.csv(csvDf, file=paste0(file.nm, ".csv"));
  
  return(1);
}
