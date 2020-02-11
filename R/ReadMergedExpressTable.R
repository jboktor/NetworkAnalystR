# read meta-dataset previously processed
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataName PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ReadMergedExpressTable
#' @export 
ReadMergedExpressTable <- function(dataName){
  current.msg <<- "";
  meta.upload <<- TRUE;
  dataSet <- .readTabData(dataName);
  common.matrix <- dataSet$data;
  meta.nms <- tolower(names(dataSet$meta.info));
  
  cls.inx <- grep("condition", meta.nms);
  if(length(cls.inx) == 0){
    current.msg <<- "No condition label found (#CLASS.condition)";
    return("F");
  }else{
    cls.inx <- cls.inx[1];
    cls.lbl <- dataSet$meta.info[[cls.inx]];
  }
  
  data.inx <- grep("dataset", meta.nms);
  if(length(data.inx) == 0){
    current.msg <<- "No dataset label found (#CLASS.dataset)";
    return("F");
  }else{
    data.inx <- data.inx[1];
    data.lbl <- dataSet$meta.info[[data.inx]];
    data.nms <- unique(as.character(data.lbl));
    
    # now create the mdata.all object
    mdata.all <- vector(mode="list", length=length(data.nms));
    names(mdata.all) <- data.nms;
    mdata.all <<- lapply(mdata.all, function(x){x=1});
  }
  
  if(length(grep("entrez.hsa", meta.nms)) > 0){
    data.org <<- "hsa"
    id.type <<- "entrez";
    shared.nms <- rownames(common.matrix);
    symbols <- doEntrez2SymbolMapping(shared.nms);
    names(symbols) <- shared.nms;
  }else if(length(grep("entrez.mmu", meta.nms)) > 0){
    data.org <<- "mmu"
    id.type <<- "entrez";
    shared.nms <- rownames(common.matrix);
    symbols <- doEntrez2SymbolMapping(shared.nms);
    names(symbols) <- shared.nms;
  }else{
    symbols <- NULL;
    inmex.org <<- "NA"
    id.type <<- 'NA';
  }
  data.org <<- unlist(strsplit( meta.nms[2], "[.]"))[3]
  inmex.meta.orig <- list(data=common.matrix,
                          id.type = id.type,
                          gene.symbls = symbols,
                          cls.lbl=factor(cls.lbl),
                          data.lbl=data.lbl);
  saveRDS(inmex.meta.orig, "inmex.meta.orig.rds");
  if(length(levels(as.factor(data.lbl))) == 1){
    return(2);
  }else{
    return(1);
  }    
}
