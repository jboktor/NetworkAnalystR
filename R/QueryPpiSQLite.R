#######################################
### Utility Methods not for public call
########################################
# note, last two par only for STRING database
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param table.nm PARAM_DESCRIPTION
#' @param q.vec PARAM_DESCRIPTION
#' @param requireExp PARAM_DESCRIPTION
#' @param min.score PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname QueryPpiSQLite
#' @export 
QueryPpiSQLite <- function(table.nm, q.vec, requireExp, min.score){
  require('RSQLite')
  ppi.db <- dbConnect(SQLite(), paste(sqlite.path, "ppi.sqlite", sep="")); 
  query <- paste(shQuote(q.vec),collapse=",");
  
  if(grepl("string$", table.nm)){
    if(requireExp){
      statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND combined_score >=", min.score, " AND experimental > 0", sep="");
    }else{
      statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND combined_score >=", min.score, sep="");
    }
  }else{
    statement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))", sep="");
  }
  ppi.res <- .query.sqlite(ppi.db, statement);
  
  # remove dupliated edges
  ppi.res <- ppi.res[!duplicated(ppi.res$row_id),]
  return(ppi.res);  
}
