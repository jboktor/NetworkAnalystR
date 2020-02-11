#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dat PARAM_DESCRIPTION
#' @param imgNm PARAM_DESCRIPTION
#' @param dpi PARAM_DESCRIPTION, Default: 72
#' @param format PARAM_DESCRIPTION, Default: 'png'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname qc.boxplot
#' @export 
qc.boxplot <- function(dat, imgNm, dpi=72, format="png"){
  dpi = as.numeric(dpi)
  library('ggplot2')
  library('lattice');
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  subgene <- 10000;
  if (nrow(dat)>subgene) {
    set.seed(28051968);
    sg  <- sample(nrow(dat), subgene);
    Mss <- dat[sg,,drop=FALSE];
  } else {
    Mss <- dat;
  }
  
  subsmpl=100;
  if (ncol(Mss)>subsmpl) {
    set.seed(28051968);
    ss  <- sample(ncol(Mss), subsmpl)
    Mss <- Mss[,ss,drop=FALSE]
  } else {
    Mss <- Mss
  }
  
  sample_id <- rep(seq_len(ncol(Mss)), each = nrow(Mss));
  values  <- as.numeric(Mss)
  
  df = cbind(values, sample_id)
  
  df = data.frame(df)
  df$sample_id = factor(df$sample_id)
  xlower = unname(quantile(df$values, probs = c(0.01, 0.99), na.rm=TRUE)[1])
  xupper = unname(quantile(df$values, probs = c(0.01, 0.99), na.rm=TRUE)[2])
  height = length(unique(df$sample_id)) *20;
  if(height<450){
    height = 450
  }
  bp = ggplot(df, aes(sample_id, values)) +
    ylab("Values") + xlab("Samples") + scale_x_discrete(labels=colnames(dataSet$data.norm)) + ylim(xlower, xupper) + stat_boxplot(geom = "errorbar", color="black")+ geom_boxplot(outlier.size=0.5, outlier.alpha=0.4)
  bp = bp + coord_flip();
  
  Cairo(file=imgNm, width=600*dpi/72, height=height*dpi/72, unit="px",dpi=dpi, type=format, bg="white");
  print(bp);
  dev.off();
}
