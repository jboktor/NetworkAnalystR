# private method for all sqlite queries
.query.sqlite <- function(db.con, statement, offline=TRUE){
  rs <- DBI::dbSendQuery(db.con, statement);
  res <- DBI::fetch(rs, n=-1); # get all records
  DBI::dbClearResult(rs);
  if(offline){
    DBI::dbDisconnect(db.con);
  }
  cleanMem();
  return(res);
}
