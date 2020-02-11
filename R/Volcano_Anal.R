#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param paired PARAM_DESCRIPTION, Default: FALSE
#' @param fcthresh PARAM_DESCRIPTION
#' @param threshp PARAM_DESCRIPTION
#' @param analType PARAM_DESCRIPTION
#' @param inx PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Volcano.Anal
#' @export 
Volcano.Anal <- function(paired=FALSE, fcthresh, threshp, analType, inx){
  inx = as.numeric(inx)
  print("Prepare volcano anal");
  if(anal.type == "metadata"){
    if(dataSet$name != selDataNm){
      dataSet <- readRDS(selDataNm);
    }
    data <- as.matrix(inmex.ind[selDataNm][[1]])
    p.value <- data[, "Pval"]
    fcthresh = 0;
    
  }else{
    data <- as.matrix(dataSet$resTable);
    
    if(dataSet$type == "array"){
      p.value <- data[, "adj.P.Val"];
    } else {
      p.value <- data[, "adj.P.Val"];
    }
  }
  fcthreshu <<- fcthresh
  
  if (analType == "qPCR"){
    inx.p <- p.value < 1;
  } else {
    inx.p <- p.value <= threshp;
  }
  zero.inx <- p.value == 0;
  if(sum(zero.inx)>0){
    p.value[zero.inx] <- min(p.value[!zero.inx])/10;
  }
  p.log <- -log10(p.value);
  
  if (dataSet$annotated){ # annotated to entrez
    anot.id <- rownames(data);
    gene.anot <- doEntrezIDAnot(anot.id);
  }else{
    anot.id <- rownames(data);
    gene.anot <- data.frame(gene_id=anot.id, symbol=anot.id, stringsAsFactors=FALSE)
    init.lib <<- "NA"
  }
  
  #gene symbol to be used for boxplot   
  
  # create a named matrix of sig vars for display
  fc.log <- data[, inx];
  hit.maxPos <- (which(fc.log> 10) )
  hit.maxNeg <- (which(fc.log< -10) )
  fc.log[hit.maxPos] = 10;
  fc.log[hit.maxNeg] = 10;
  #fc.all <- res$fc.all;
  
  if(fcthresh != 0){
    inx.up = fc.log > fcthresh & p.value < threshp;
    inx.down = fc.log < -fcthresh & p.value < threshp;
  }else{
    inx.up = fc.log > 0 & p.value < threshp;
    inx.down = fc.log < 0 & p.value < threshp;
  }
  
  # create named sig table for display
  inx.imp <- (inx.up | inx.down) & inx.p;
  sig.var <- cbind(fc.log[inx.imp,drop=F], p.value[inx.imp,drop=F], p.log[inx.imp,drop=F]);
  colnames(sig.var) <- c("log2(FC)", "p.value", "-log10(p)");
  # first order by log(p), then by log(FC)
  ord.inx <- order(sig.var[,3], abs(sig.var[,1]), decreasing=T);
  sig.var <- sig.var[ord.inx,, drop=F];
  
  sig.var <- signif (sig.var, 5);
  sig.var1 <- sig.var;
  sig.var1 = cbind(rownames(sig.var), sig.var);
  colnames(sig.var1) <- c("name", "log2(FC)", "p.value", "-log10(p)");
  
  ###########################
  ## for Volcano data
  ##########################
  
  if(init.lib != "NA"){
    PerformVolcanoEnrichment("abc", init.lib, "null", "all", inx)
  }
  
  fileName <- "volcano.csv";
  jsonNm <- "volcano.json";
  require(RJSONIO);
  json.obj <- toJSON(sig.var1);
  sink(jsonNm);
  cat(json.obj);
  sink();
  write.csv(signif (sig.var,5),file=fileName);
  colnames(gene.anot)[1] = "anot.id"
  volcano <- list (
    raw.threshx = fcthresh,
    raw.threshy = threshp,
    paired = paired,
    thresh.y = -log10(threshp),
    fc.symb =rownames(data),
    fc.log = fc.log,
    fc.log.uniq = jitter(fc.log),
    inx.up = inx.up,
    inx.down = inx.down,
    p.log = p.log,
    inx.p = inx.p,
    sig.mat = sig.var,
    conv = gene.anot
  );
  
  require(RJSONIO);
  json.obj <- toJSON(volcano);
  sink("volcano2.json");
  cat(json.obj);
  sink();
  
  require(RJSONIO);
  
  if(init.lib == "NA"){
    enr.mat = "NA"
  }
  write.csv(enr.mat, file="enrichment_result.csv", row.names=T);
  sink("enrichment_result.json");
  cat(json.obj);
  sink();
  
}
