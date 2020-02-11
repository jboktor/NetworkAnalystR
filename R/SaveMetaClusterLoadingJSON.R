#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fileName PARAM_DESCRIPTION
#' @param clustOpt PARAM_DESCRIPTION
#' @param nb PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SaveMetaClusterLoadingJSON
#' @export 
SaveMetaClusterLoadingJSON <- function(fileName, clustOpt, nb){
  
  inmex.meta <- readRDS("inmex_meta.rds");
  datanm.vec <- names(mdata.all)[mdata.all==1];
  nb = as.numeric(nb)
  dat.inx <- inmex.meta$data.lbl %in% datanm.vec;
  dat <- inmex.meta$data[, dat.inx, drop=F]; 
  
  # need to deal with missing values 
  dat <- na.omit(dat);
  variances = apply(dat,1, function(x){var(x)})
  df = data.frame(var = variances, inx = seq.int(1,length(variances)))
  df = df[order(-df$var),]
  inx = df$inx[c(1:nb)]
  dat = dat[inx,];
  
  pca3d <- list();
  
  pca <- prcomp(t(dat), center=T, scale=T);    
  imp.pca<-summary(pca)$importance;
  pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
  coords <- data.frame(t(signif(pca$rotation[,1:3], 5)));
  
  colnames(coords) <- NULL; 
  pca3d$score$xyz <- coords;
  pca3d$score$name <- doEntrez2SymbolMapping(rownames(pca$rotation));
  pca3d$score$entrez <- rownames(pca$rotation);
  
  loadEntrez <<- pca3d$score$entrez
  mypos <- t(coords);
  colnames(mypos) <- paste("Dim", 1:3, sep="");
  
  coords <- data.frame(mypos);
  write.csv(coords, file="networkanalyst_loadings_3d_pos.csv");
  
  require(RJSONIO);
  json.mat <- toJSON(pca3d, .na='null');
  sink(fileName);
  cat(json.mat);
  sink();
  current.msg <<- "Annotated data is now ready for 3D visualization!";
  return(1);
}
