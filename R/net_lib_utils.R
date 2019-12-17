##################################################
## R script for NetworkAnalyst
## Description: Functions to load various libraries for functional enrichment analysis during network visualization
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################

# table name is org code, id.type is column name
QueryMirSQLite <- function(org, id.type, q.vec){
  require('RSQLite');
  mir.db <- dbConnect(SQLite(), paste(sqlite.path, "mir.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", org, " WHERE ",id.type," IN (",query,")", sep="");
  mirtable <- dbSendQuery(mir.db, statement);
  mir.dic <- fetch(mirtable, n=-1); # get all records
  dbDisconnect(mir.db);
  return(mir.dic)}

# table name is org code, id.type is column name
QueryDrugSQLite <- function(q.vec){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, "drug.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM human WHERE upid IN (",query,")", sep="");
  drugtable <- dbSendQuery(drug.db, statement);
  drug.dic <- fetch(drugtable, n=-1); # get all records
  dbDisconnect(drug.db);
  return(drug.dic)}

QueryDiseaseSQLite <- function(q.vec){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, "disease.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM human WHERE entrez IN (",query,")", sep="");
  drugtable <- dbSendQuery(drug.db, statement);
  drug.dic <- fetch(drugtable, n=-1); # get all records
  dbDisconnect(drug.db);
  return(drug.dic)}

QueryTfmirSQLite <- function(q.vec){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, "tfmir.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", data.org, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))" , sep="");
  drugtable <- dbSendQuery(drug.db, statement);
  drug.dic <- fetch(drugtable, n=-1); # get all records
  dbDisconnect(drug.db);
  return(drug.dic)}

QueryDiffNetSQLite <- function(q.vec){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, "tissuePPI.sqlite", sep=""));
  table.nm = diffNetName;
  query <- paste (shQuote(q.vec),collapse=",");
  topPct = 1-diffPct;
  botstatement <- paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND rank <=", diffPct ,sep="");
  paste("SELECT * FROM ",table.nm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))  AND rank >=", topPct ,sep="") ;
  drugtablebot <- dbSendQuery(drug.db, botstatement);
  drug.dic1 <- fetch(drugtablebot, n=-1); # get all records
  dbClearResult(drugtablebot);
  drugtabletop <- dbSendQuery(drug.db, botstatement);
  drug.dic2 <- fetch(drugtabletop, n=-1); # get all records
  drug.dic = rbind(drug.dic1, drug.dic2);
  dbDisconnect(drug.db);
  return(drug.dic)}

QueryCellCoexSQLite <- function(q.vec){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, data.org,"_immune.sqlite", sep=""));
  tblNm = paste0(data.org,"_",cellCoexNumber);
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", tblNm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))" , sep="");
  drugtable <- dbSendQuery(drug.db, statement);
  drug.dic <- fetch(drugtable, n=-1); # get all records
  dbDisconnect(drug.db);
  return(drug.dic)}

QueryTissueCoexSQLite <- function(q.vec){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, "tissueCoex.sqlite", sep=""));
  tblNm = paste0(data.org,"_",cellCoexNumber);
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", tblNm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))" , sep="");
  drugtable <- dbSendQuery(drug.db, statement);
  drug.dic <- fetch(drugtable, n=-1); # get all records
  dbDisconnect(drug.db);
  return(drug.dic)}

QueryChemSQLite<- function(org, q.vec){
  require('RSQLite');
  chem.db <- dbConnect(SQLite(), paste(sqlite.path, "chem.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", org, " WHERE entrez IN (",query,")", sep="");
  chemtable <- dbSendQuery(chem.db, statement);
  chem.dic <- fetch(chemtable, n=-1); # get all records
  dbDisconnect(chem.db);
  return(chem.dic)}

QueryTFSQLite<- function(table.nm, q.vec){
  require('RSQLite');
  chem.db <- dbConnect(SQLite(), paste(sqlite.path, "tfac.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", table.nm, " WHERE entrez IN (",query,")", sep="");
  chemtable <- dbSendQuery(chem.db, statement);
  chem.dic <- fetch(chemtable, n=-1); # get all records
  dbDisconnect(chem.db);
  return(chem.dic)}

doPpiIDMapping <- function(q.vec){
  if(data.org == "ath"){
    db.path <- paste(lib.path, data.org, "/tair.rds", sep="");
  }else if(data.org == "sce"){
    if(net.type == "string"){ # only for yeast
      db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
    }else{
      db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep="");
    }
  }else if(net.type %in% c("innate", "irefinx", "rolland")){
    db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep="");
  }else{
    if(data.org %in% c("bsu", "tbr", "cel", "dme", "eco", "pfa","pae")){
      db.path <- paste(lib.path, data.org, "/entrez_string.rds", sep="");
    }else if (data.org %in% c("mmu","hsa") && net.type == "string"){
      db.path <- paste(lib.path, data.org, "/string.rds", sep="");
    }else if(data.org == "mtb"){
      db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
    }else{
      db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
    }
  }
  db.map <-  readRDS(db.path);
  hit.inx <- match(q.vec, db.map[, "gene_id"]);
  ppi.mat <- db.map[hit.inx, ];

  # fix the factor col related to library issue
  i <- sapply(ppi.mat, is.factor)
  ppi.mat[i] <- lapply(ppi.mat[i], as.character)
  if(data.org %in% c("pae", "mtb")){
    ppi.mat = ppi.mat[,c(2,1)];
    colnames(ppi.mat) = c("gene_id", "accession");
  }
  return(ppi.mat)}

doUniprot2EntrezMapping <- function(uniprot.vec){
  db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep="");
  db.map <-  readRDS(db.path);
  hit.inx <- match(uniprot.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs)}

doEntrez2UniprotMapping <- function(entrez.vec){
  db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep="");
  db.map <-  readRDS(db.path);
  hit.inx <- match(entrez.vec, db.map[, "gene_id"]);
  unips <- db.map[hit.inx, "accession"];
  mode(unips) <- "character";
  return(unips)}

doString2EntrezMapping <- function(string.vec){
  db.path <- paste(lib.path, data.org, "/entrez_string.rds", sep="");
  db.map <-  readRDS(db.path);
  hit.inx <- match(string.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs)}

doEmblGene2EntrezMapping <- function(emblgene.vec){
  db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
  db.map <-  readRDS(db.path);
  hit.inx <- match(emblgene.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs)}

doEmblProtein2EntrezMapping <- function(emblprotein.vec){
  db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
  db.map <-  readRDS(db.path);
  hit.inx <- match(emblprotein.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs)}
