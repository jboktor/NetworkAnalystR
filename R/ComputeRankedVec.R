#For GSEA of AnalOverview page
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param inx PARAM_DESCRIPTION, Default: 1
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ComputeRankedVec
#' @export 
ComputeRankedVec <- function(data, inx = 1){
  opt = rankOptGlobal;
  cls = data$cls
  if(anal.type == "metadata"){
    if(opt %in% c("mwt", "s2n", "wcx", "stu")){
      matr <- as.matrix(data$data)
    }else{
      matr <- as.matrix(readRDS("meta.resTable.rds"));
    }
  }else{
    if(opt %in% c("mwt", "s2n", "wcx", "stu")){
      matr = as.matrix(data$data.norm)
    }else{
      matr = as.matrix(data$resTable);
    }
  }
  if(opt == "mwt"){
    res = CalculateMWT(matr, cls)
    rankedVec = res$MWT
    names(rankedVec) = rownames(matr)
  }else if(opt == "stu"){
    inx1 <- which(data$cls==levels(data$cls)[1]);
    inx2 <- which(data$cls==levels(data$cls)[2]);
    res <- apply(matr, 1, function(x) {
      tmp <- try(t.test(x[inx1], x[inx2], paired = FALSE, var.equal = TRUE));
      if(class(tmp) == "try-error") {
        return(c(NA, NA));
      }else{
        return(c(tmp$statistic, tmp$p.value));
      }
    })
    res = t(res)
    rankedVec = res[,1]
    posInx = sign(rankedVec) == 1
    rankedVec[posInx] = 1000-rankedVec[posInx]
    names(rankedVec) = rownames(matr)
    
  }else if(opt == "wcx"){
    inx1 <- which(data$cls==levels(data$cls)[1]);
    inx2 <- which(data$cls==levels(data$cls)[2]);
    res <- apply(matr, 1, function(x) {
      tmp <- try(wilcox.test(x[inx1], x[inx2], paired = FALSE));
      if(class(tmp) == "try-error") {
        return(c(NA, NA));
      }else{
        return(c(tmp$statistic, tmp$p.value));
      }
    })
    res = t(res)
    rankedVec = res[,2]
    names(rankedVec) = rownames(matr)
    
  }else if (opt == "s2n"){
    res = CalculateS2N(matr, as.numeric(cls)-1)
    rankedVec = res;
  }else if(opt == "pval"){
    m = length(colnames(matr))
    if (dataSet$de.method=="limma"){
      rankedVec = as.vector(matr[,"t"]);
    } else if (dataSet$de.method=="deseq2"){
      rankedVec = as.vector(matr[,"stat"]);
    } else {
      rankedVec = as.vector(matr[,"LR"]);
    }
    names(rankedVec) = rownames(matr);
  }else{
    rankedVec = as.vector(matr[,inx]);
    names(rankedVec) = rownames(matr);
  }
  rankedVec = sort(rankedVec)
  rankedVec = rankedVec[unique(names(rankedVec))]
  rankedVec <<- rankedVec
  return(rankedVec)
}
