#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cmpdNm PARAM_DESCRIPTION
#' @param dpi PARAM_DESCRIPTION, Default: 72
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PlotMetaPhm
#' @export 
PlotMetaPhm <-function(cmpdNm, dpi=72){
  allmat = readRDS("allmat.rds");
  
  fileNm = paste("Path_", cmpdNm, ".png", sep="");
  current.geneset <- readRDS("current_geneset.rds")
  ids=current.geneset[[cmpdNm]];
  subset = allmat[which(rownames(allmat) %in% ids),]
  if(length(subset)<1){
    subset = allmat[which(rownames(allmat) %in% ids),]
  }
  library(RColorBrewer);
  library(pheatmap)
  dims = dim(subset)
  rnms = rownames(subset)
  cnms = colnames(subset) 
  m <- mapply(subset, FUN=as.numeric)
  subset <- matrix(data=m, ncol=dims[2], nrow=dims[1])
  rownames(subset) = rnms
  colnames(subset) = cnms
  if(inmex.method %in% c("effectsize", "merge")){
    subset = subset[complete.cases(subset), ];
  }
  
  subset[is.na(subset)] <- 0
  my_palette <- colorRampPalette(c("green", "black", "red"))(n = 30)
  my_palette = c("#d3d3d3", my_palette)
  bk2 = unique(c(seq(0,0.9999, length=2), 1, seq(2,30, length=29)));
  inmex.meta <- readRDS("inmex_meta.rds");
  ann = data.frame(Class=inmex.meta$cls.lbl, Dataset=inmex.meta$data.lbl)
  rownames(ann) = colnames(subset)
  Cairo(file=fileNm, width=800, height=700, type="png", bg="white",unit="px",dpi=72)
  hm = pheatmap(subset, color = my_palette, breaks = bk2, annotation_col=ann, show_rownames = FALSE, main = cmpdNm, border_color=NA)
  print(hm)
  dev.off()
}
