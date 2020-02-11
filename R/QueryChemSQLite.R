#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param org PARAM_DESCRIPTION
#' @param q.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname QueryChemSQLite
#' @export 
QueryChemSQLite<- function(org, q.vec){
  require('RSQLite');
  chem.db <- dbConnect(SQLite(), paste(sqlite.path, "chem.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", org, " WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(chem.db, statement));
}
