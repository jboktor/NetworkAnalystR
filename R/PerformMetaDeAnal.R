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
#' @rdname PerformMetaDeAnal
#' @export 
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
