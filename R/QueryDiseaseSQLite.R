#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param q.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname QueryDiseaseSQLite
#' @export 
QueryDiseaseSQLite <- function(q.vec){
  require('RSQLite');
  dis.db <- dbConnect(SQLite(), paste(sqlite.path, "disease.sqlite", sep="")); 
  query <- paste(shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM human WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(dis.db, statement));
}
