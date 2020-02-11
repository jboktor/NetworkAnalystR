# for gene set-level meta-analysis
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @param netNm PARAM_DESCRIPTION
#' @param BHth PARAM_DESCRIPTION
#' @param mType PARAM_DESCRIPTION
#' @param curr.geneset PARAM_DESCRIPTION
#' @param lib PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[RJSONIO]{toJSON}}
#' @rdname SetupMetaGSEAStats
#' @export 
#' @importFrom RJSONIO toJSON
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
