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
#' @rdname QueryDiffNetSQLite
#' @export 
QueryDiffNetSQLite <- function(q.vec){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, "tissuePPI.sqlite", sep="")); 
  table.nm = diffNetName;
  query <- paste (shQuote(q.vec),collapse=",");
  topPct = 1-diffPct;
  botstatement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND rank <=", diffPct ,sep="");
  topstatement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND rank >=", topPct ,sep="");
  
  drug.dic1 <- .query.sqlite(drug.db, botstatement, FALSE);# no close db connection
  drug.dic2 <- .query.sqlite(drug.db, topstatement);
  drug.dic <- rbind(drug.dic1, drug.dic2);
  return(drug.dic);
}
