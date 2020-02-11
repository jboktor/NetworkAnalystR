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
#' @rdname SaveExpressClusterLoadingJSON
#' @export 
SaveExpressClusterLoadingJSON <- function(fileName, clustOpt, nb){
  dat <- dataSet$data.norm;
  pca3d <- list();
  dat <- na.omit(dat);
  nb = as.numeric(nb)
  if(clustOpt == "pca"){
    pca <- prcomp(t(dat), center=T, scale=T);
    imp.pca<-summary(pca)$importance;
    pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
    coords <- data.frame(t(signif(pca$rotation[,1:3], 5)));
    
    colnames(coords) <- NULL; 
    pca3d$score$xyz <- coords;
    pca3d$score$name <- doEntrez2SymbolMapping(rownames(pca$rotation));
    pca3d$score$entrez <-rownames(pca$rotation);
    weights = imp.pca[2,][1:3]
    mypos <- t(coords);
    meanpos = apply(abs(mypos),1, function(x){weighted.mean(x, weights)})
    df = data.frame(pos = meanpos, inx = seq.int(1,length(meanpos)))
    df = df[order(-df$pos),]
    if(nrow(df) >2000){
      inx = df$inx[c(1:nb)]
      mypos = mypos[inx,];
      pca3d$score$xyz = coords[inx]
      pca3d$score$name = pca3d$score$name[inx]
      pca3d$score$entrez = pca3d$score$entrez[inx]
    }
  }
  
  pca3d$cls = dataSet$meta.info;
  colnames(mypos) <- paste("Dim", 1:3, sep="");
  # see if there is secondary
  loadEntrez <<- pca3d$score$entrez
  rownames(mypos) = pca3d$score$name;
  
  write.csv(mypos, file="networkanalyst_3d_load_pos.csv");
  require(RJSONIO);
  json.mat <- toJSON(pca3d, .na='null');
  sink(fileName);
  cat(json.mat);
  sink();
  current.msg <<- "Annotated data is now ready for PCA 3D visualization!";
  return(1);
}
