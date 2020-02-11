#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param table.nm PARAM_DESCRIPTION
#' @param q.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname QueryTFSQLite
#' @export 
QueryTFSQLite<- function(table.nm, q.vec){
  require('RSQLite');
  chem.db <- dbConnect(SQLite(), paste(sqlite.path, "tfac.sqlite", sep="")); 
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", table.nm, " WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(chem.db, statement));
}
