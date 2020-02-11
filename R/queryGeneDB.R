#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param table.nm PARAM_DESCRIPTION
#' @param data.org PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname queryGeneDB
#' @export 
queryGeneDB <- function(table.nm, data.org){
  require('RSQLite')
  
  conv.db <- dbConnect(SQLite(), paste(genesdb.path, data.org, "_genes.sqlite", sep="")); 
  db.map <- dbReadTable(conv.db, table.nm);
  dbDisconnect(conv.db); cleanMem();
  
  return(db.map)
}
