# utils to remove from
# within, leading and trailing spaces
# remove /
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param query PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ClearStrings
#' @export 
ClearStrings<-function(query){
  # remove leading and trailing space
  query<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", query, perl=TRUE);
  
  # kill multiple white space
  query <- gsub(" +",".",query);
  query <- gsub("/", ".", query);
  query <- gsub("-", ".", query);
  return (query);
}
