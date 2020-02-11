# retrun the json obj
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
#' @rdname SaveMetaClusterJSON
#' @export 
SaveMetaClusterJSON <- function(fileName, clustOpt, opt){
  
  inmex.meta <- readRDS("inmex_meta.rds");
  datanm.vec <- names(mdata.all)[mdata.all==1];
  dat.inx <- inmex.meta$data.lbl %in% datanm.vec;
  dat <- inmex.meta$data[, dat.inx, drop=F]; 
  
  # need to deal with missing values 
  dat <- na.omit(dat);
  
  pca3d <- list();
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
  }else{
    require('Rtsne');
    ndat <- as.matrix(t(dat));
    max.perx <- floor((nrow(ndat)-1)/3);
    if(max.perx > 30){
      max.perx <- 30;
    }
    res <- Rtsne(ndat, dims = 3, perplexity=max.perx);
    pca3d$score$axis <- paste("t-SNE dim ", 1:3, sep="");
    coords <- data.frame(t(signif(res$Y, 5)));
  }
  
  colnames(coords) <- NULL; 
  pca3d$score$xyz <- coords;
  pca3d$score$name <- colnames(dat);
  
  facA <- as.character(inmex.meta$cls.lbl[dat.inx]);
  if(all.numeric(facA)){
    facA <- paste("Group", facA);
  }
  pca3d$score$facA <- facA;
  
  facB <-  as.character(inmex.meta$data.lbl[dat.inx]);
  if(all.numeric(facB)){
    facB <- paste("Group", facB);
  }
  pca3d$score$facB <- facB;
  
  # now set color for each group
  cols <- unique(GetColorSchema(facB));
  rgbcols <- col2rgb(cols);
  cols <- apply(rgbcols, 2, function(x){paste("rgba(", paste(x, collapse=","), ",1)", sep="")});
  pca3d$score$colors <- cols;
  
  # add shape sphere, triangles, square, pentagon (first two)
  pca3d$score$shapes <- c("sphere", "triangle");
  
  mypos <- t(coords);
  colnames(mypos) <- paste("Dim", 1:3, sep="");
  coords <- data.frame(Class=facA, Data=facB, mypos);
  write.csv(coords, file="networkanalyst_3d_pos.csv");
  
  require(RJSONIO);
  json.mat <- toJSON(pca3d, .na='null');
  sink(fileName);
  cat(json.mat);
  sink();
  current.msg <<- "Annotated data is now ready for 3D visualization!";
  return(1);
}
