#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fileName PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[data.table]{fread}}
#' @rdname ReadListData
#' @export 
#' @importFrom data.table fread
ReadListData <- function(fileName) {
  dat1 <- data.table::fread(fileName, header=FALSE, check.names=FALSE, data.table=FALSE);
  dataSet$name = fileName
  rowNms = dat1[,1]
  if(length(dat1) == 1){
    dat1[,1] = 0
  }else{
    dat1[,1] = dat1[,2]
    dat1 = dat1[,-2];
  }
  dataSet$prot.mat = as.matrix(dat1)
  rownames(dataSet$prot.mat) = rowNms;
  saveRDS(dataSet, file=fileName); # keep original copy, not in mem
  return(1)
}
