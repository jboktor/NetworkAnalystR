#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetDensityPlot
#' @export 
GetDensityPlot <- function() {
  data = list();
  dat = dataSet$data.norm;
  for(i in 1: ncol(dat)){
    data[[i]] = dat[,i]
    names(data)[i] = colnames(dat)[i];
  }
  densityList = list();
  for(i in 1: length(data)){
    d = density(data[[i]])
    df = data.frame(d$x, d$y)
    colnames(df) = c("x","y");
    densityList[[i]] = df
  }
  names(densityList) = colnames(dataSet$data.norm);
  jsonNm <- "density.json";
  lst = list()
  
  lst=list(
    density= densityList,
    class= dataSet$cls
  )
  require(RJSONIO);
  json.obj <- toJSON(lst);
  sink(jsonNm);
  cat(json.obj);
}
