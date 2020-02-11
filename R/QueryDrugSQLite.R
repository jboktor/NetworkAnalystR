# table name is org code, id.type is column name
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
#' @rdname QueryDrugSQLite
#' @export 
QueryDrugSQLite <- function(q.vec){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, "drug.sqlite", sep="")); 
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM human WHERE upid IN (",query,")", sep="");
  return(.query.sqlite(drug.db, statement));
}
