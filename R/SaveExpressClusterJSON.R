# single expression data
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fileName PARAM_DESCRIPTION
#' @param clustOpt PARAM_DESCRIPTION
#' @param opt PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SaveExpressClusterJSON
#' @export 
SaveExpressClusterJSON <- function(fileName, clustOpt, opt){
  dat <- dataSet$data.norm;
  pca3d <- list();
  dat <- na.omit(dat);
  
  if(clustOpt == "pca"){
    if(opt == "all"){
      pca <- prcomp(t(dat), center=T, scale=T);
    }else{
      dat = dat[which(rownames(dat) %in% loadEntrez),]
      pca <- prcomp(t(dat), center=T, scale=T);
    }
    imp.pca<-summary(pca)$importance;
    pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
    coords <- data.frame(t(signif(pca$x[,1:3], 5)));
  }else{ # tsne
    require('Rtsne');
    dat <- as.matrix(t(dat));
    max.perx <- floor((nrow(dat)-1)/3);
    if(max.perx > 30){
      max.perx <- 30;
    }
    res <- Rtsne(dat, dims = 3, perplexity=max.perx);
    pca3d$score$axis <- paste("t-SNE dim ", 1:3, sep="");
    coords <- data.frame(t(signif(res$Y, 5)));
  }
  
  colnames(coords) <- NULL; 
  pca3d$score$xyz <- coords;
  pca3d$score$name <- colnames(dataSet$data.norm);
  
  facA <- as.character(dataSet$fst.cls);
  if(all.numeric(facA)){
    facA <- paste("Group", facA);
  }
  pca3d$score$facA <- facA;
  
  mypos <- t(coords);
  colnames(mypos) <- paste("Dim", 1:3, sep="");
  # see if there is secondary
  if(length(dataSet$sec.cls) > 1){
    facB <- as.character(dataSet$sec.cls);
    if(all.numeric(facB)){
      facB <- paste("Group", facB);
    }
    pca3d$score$facB <- facB;
    
    # set shape based on the first group
    pca3d$score$shapes <- c("sphere", "triangle");
    
    # now set color based on 2nd group
    cols <- unique(GetColorSchema(dataSet$sec.cls));
    rgbcols <- col2rgb(cols);
    cols <- apply(rgbcols, 2, function(x){paste("rgb(", paste(x, collapse=","), ")", sep="")});
    pca3d$score$colors <- cols;
    
    mypos <- data.frame(factorA=facA, factorB=facB, mypos);
  }else{
    # now set color based on first group
    cols <- unique(GetColorSchema(dataSet$fst.cls));
    rgbcols <- col2rgb(cols);
    cols <- apply(rgbcols, 2, function(x){paste("rgba(", paste(x, collapse=","), ",1)", sep="")});
    pca3d$score$colors <- cols;
    mypos <- data.frame(factorA=facA, mypos);
  }
  pca3d$cls = dataSet$meta.info;
  rownames(mypos) = colnames(dataSet$data.norm);
  
  write.csv(mypos, file="networkanalyst_3d_pos.csv");
  require(RJSONIO);
  json.mat <- toJSON(pca3d, .na='null');
  sink(fileName);
  cat(json.mat);
  sink();
  current.msg <<- "Annotated data is now ready for PCA 3D visualization!";
  return(1);
}
