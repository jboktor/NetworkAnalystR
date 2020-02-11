# note, we do both filtering and normalization
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param norm.opt PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformDataNormalization
#' @export 
PerformDataNormalization <- function(data, norm.opt){
  set.seed(1337);
  msg <- NULL;
  row.nms <- rownames(data);
  col.nms <- colnames(data);
  if(norm.opt=="log"){
    min.val <- min(data[data>0], na.rm=T)/10;
    numberOfNeg = sum(data<=0, na.rm = TRUE) + 1; 
    totalNumber = length(data)
    if((numberOfNeg/totalNumber)>0.2){
      msg <- paste(msg, "Can't perform log2 normalization, over 20% of data are negative. Try a different method or maybe the data already normalized?", collapse=" ");
      print(msg);
      norm.msg <<- current.msg <<- msg;
      return(0);
    }
    data[data<=0] <- min.val;
    data <- log2(data);
    msg <- paste(msg, "Log2 transformation.", collapse=" ");
  }else if(norm.opt=="vsn"){
    require(limma);
    data <- normalizeVSN(data);
    msg <- paste(msg, "VSN normalization.", collapse=" ");
  }else if(norm.opt=="quantile"){
    require('preprocessCore');
    data <- normalize.quantiles(data, copy=TRUE);
    msg <- paste(msg, "Quantile normalization.", collapse=" ");
  }else if(norm.opt=="combined"){
    require(limma);
    data <- normalizeVSN(data);
    require('preprocessCore');
    data <- normalize.quantiles(data, copy=TRUE);
    msg <- paste(msg, "VSN followed by quantile normalization.", collapse=" ");
  }else if(norm.opt=="logcount"){ # for count data, do it in DE analysis, as it is dependent on design matrix
    require(edgeR);
    nf <- calcNormFactors(data);
    y <- voom(data,plot=F,lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    msg <- paste(msg, "Limma based on log2-counts per million transformation.", collapse=" ");
  }else{
    # should do best guess for count data for plotting and filtering
    if(dataSet$type == "count"){
      if(sum(data > 100) > 100){ # now we think it is raw counts
        require(edgeR);
        nf <- calcNormFactors(data);
        y <- voom(data,plot=F,lib.size=colSums(data)*nf);
        data <- y$E; # copy per million
      }
    }
    msg <- paste(msg, "No log normalization was performed.", collapse=" ");
    print(msg);
  }
  norm.msg <<- msg;
  rownames(data) <- row.nms;
  colnames(data) <- col.nms;
  return(data);
}
