#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param norm.opt PARAM_DESCRIPTION
#' @param var.thresh PARAM_DESCRIPTION
#' @param count.thresh PARAM_DESCRIPTION
#' @param abundance PARAM_DESCRIPTION
#' @param filterUnmapped PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformExpressNormalization
#' @export 
PerformExpressNormalization <- function(norm.opt, var.thresh, count.thresh, abundance, filterUnmapped){
  
  print("normalizing ....");
  msg <- "Only features with annotations are kept for further analysis.";
  
  if(filterUnmapped == "false"){
    # need to update those with annotations
    data1 <- readRDS("data.proc.rds");
    anot.id <- readRDS("annotation.rds");;
    hit.inx <- !is.na(anot.id);
    rownames(data1)[hit.inx] <- anot.id[hit.inx];
    data1 <- RemoveDuplicates(data1, "mean", quiet=T);
    raw.data.anot <- data <- dataSet$data.anot <- data1;
  }else{
    raw.data.anot <- data <- readRDS("orig.data.anot");
  }
  
  if (dataSet$type == "count"){
    sum.counts <- apply(data, 1, sum, na.rm=TRUE);
    rm.inx <- sum.counts < count.thresh;
    data <- data[!rm.inx,];
    msg <- paste(msg, "Filtered ", sum(rm.inx), " genes with low counts.", collapse=" ");
  }else{
    avg.signal <- apply(data, 1, mean, na.rm=TRUE)
    p05 <- quantile(avg.signal, 0.05)
    all.rows = nrow(data)
    rm.inx = avg.signal < p05
    data <- data[!rm.inx,]
    msg <- paste(msg, "Filtered ", sum(rm.inx), " genes with low relative abundance (average expression signal).", collapse=" ");
  }
  
  data <- PerformDataNormalization(data, norm.opt);
  if(length(data)==1 && data == 0){
    return(0);
  }
  msg <- paste(norm.msg, msg);
  
  filter.val <- apply(data, 1, IQR, na.rm=T);
  nm <- "Interquantile Range";
  rk <- rank(-filter.val, ties.method='random');
  kp.pct <- (100 - var.thresh)/100;
  
  remain <- rk < nrow(data)*kp.pct;
  data <- data[remain,];
  msg <- paste(msg, paste("Filtered ", sum(!remain), " low variance genes based on IQR"), collapse=" ");
  
  dataSet$data.anot <- raw.data.anot[remain,]
  dataSet$data.norm <- data;
  
  # save normalized data for download user option
  write.csv(dataSet$data.norm, file="data_normalized.csv");
  
  current.msg <<- msg;
  dataSet <<- dataSet;
  
  processedObj <- list();#for omicsnet
  processedObj$name <- "rna_b_omicsanalyst.json"
  processedObj$type <- "rna.b"
  processedObj$data.proc <- dataSet$data.norm
  processedObj$feature.nms <- rownames(processedObj$data.proc)
  processedObj$sample.nms <- colnames(processedObj$data.proc)
  meta = dataSet$meta.info
  rownames(meta) <-  colnames(processedObj$data.proc)
  processedObj$meta <- meta
  library(RJSONIO)
  sink(processedObj$name);
  cat(toJSON(processedObj));
  sink();
  
  return(1);
}
