# users can select one or more data for analysis
# note, we use 1 to indicate this is selected
# and by default is all selected. 
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
#' @rdname SelectData
#' @export 
SelectData <- function(){
  if(!exists('nm.vec')){
    current.msg <<-"No dataset is selected for analysis!";
    print(current.msg);
    return(0);
  }
  
  all.nms <- names(mdata.all);
  for(nm in all.nms){
    if(nm %in% nm.vec){
      mdata.all[[nm]] <<- 1;
    }else{
      mdata.all[[nm]] <<- 0;
    }
  }
  if(anal.type == "metadata"){
    if("meta_dat" %in% nm.vec){
      meta.selected <<- TRUE;
    }else{
      meta.selected <<- FALSE;
    }
  }
  rm('nm.vec', envir = .GlobalEnv);
  return(1);
}
