##################################################
## R script for NetworkAnalyst
## Description: Functions to load various libraries for functional enrichment analysis during network visualization
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################
# table name is org code, id.type is column name
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param org PARAM_DESCRIPTION
#' @param id.type PARAM_DESCRIPTION
#' @param q.vec PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname QueryMirSQLite
#' @export 
QueryMirSQLite <- function(org, id.type, q.vec){
  require('RSQLite');
  mir.db <- dbConnect(SQLite(), paste(sqlite.path, "mir.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", org, " WHERE ",id.type," IN (",query,")", sep="");
  return(.query.sqlite(mir.db, statement));
}
