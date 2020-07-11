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
  mir.db <- dbConnect(SQLite(), paste(sqlite.path, "mir2gene.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", org, " WHERE ",id.type," IN (",query,")", sep="");
  return(.query.sqlite(mir.db, statement));
}

# table name is org code, id.type is column name
QueryDrugSQLite <- function(q.vec){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, "drug.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM human WHERE upid IN (",query,")", sep="");
  return(.query.sqlite(drug.db, statement));
}

QueryDiseaseSQLite <- function(q.vec){
  require('RSQLite');
  dis.db <- dbConnect(SQLite(), paste(sqlite.path, "disease.sqlite", sep=""));
  query <- paste(shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM human WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(dis.db, statement));
}

QueryTfmirSQLite <- function(q.vec){
  require('RSQLite');
  tf.db <- dbConnect(SQLite(), paste(sqlite.path, "tfmirgene.sqlite", sep=""));
  query <- paste(shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", data.org, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))" , sep="");
  return(.query.sqlite(tf.db, statement));
}

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

QueryCellCoexSQLite <- function(q.vec){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, data.org,"_immune.sqlite", sep=""));
  tblNm = paste0(data.org,"_",cellCoexNumber);
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", tblNm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))" , sep="");
  return(.query.sqlite(drug.db,statement));
}

QueryTissueCoexSQLite <- function(q.vec){
  require('RSQLite');
  drug.db <- dbConnect(SQLite(), paste(sqlite.path, "tissueCoex.sqlite", sep=""));
  tblNm = paste0(data.org,"_",cellCoexNumber);
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", tblNm, " WHERE ((id1 IN (", query, ")) OR (id2 IN (", query, ")))" , sep="");
  return(.query.sqlite(drug.db, statement));
}

QueryChemSQLite<- function(org, q.vec){
  require('RSQLite');
  chem.db <- dbConnect(SQLite(), paste(sqlite.path, "chem.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", org, " WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(chem.db, statement));
}

QueryTFSQLite<- function(table.nm, q.vec){
  require('RSQLite');
  chem.db <- dbConnect(SQLite(), paste(sqlite.path, "tf2gene.sqlite", sep=""));
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", table.nm, " WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(chem.db, statement));
}

doPpiIDMapping <- function(q.vec, table.nm="entrez_swissprot"){
  if(data.org == "ath"){
    db.map <-  queryGeneDB("tair", data.org);
  }else if(data.org == "sce"){
    if(net.type == "string"){ # only for yeast
      db.map <-  queryGeneDB("entrez_embl_gene", data.org);
    }else{
      db.map <-  queryGeneDB("entrez_uniprot", data.org);
    }
  }else if(net.type == "innate"){
    db.map <-  queryGeneDB(table.nm, data.org);
  }else if(net.type %in% c("irefinx", "rolland")){
    db.map <-  queryGeneDB("entrez_uniprot", data.org);
  }else{
    if(data.org %in% c("bsu", "tbr", "cel", "dme", "eco", "pfa","pae")){
      db.map <-  queryGeneDB("entrez_string", data.org);
    }else if (data.org %in% c("mmu","hsa") && net.type %in% c("string", "huri")){
      db.map <-  queryGeneDB("entrez", data.org);
    }else if(data.org == "mtb"){
      db.map <-  queryGeneDB("entrez", data.org);
    }else{
      db.map <-  queryGeneDB("entrez_embl_protein", data.org);
    }
  }
  hit.inx <- match(q.vec, db.map[, "gene_id"]);
  ppi.mat <- db.map[hit.inx, ];
  if(data.org %in% c("pae", "mtb") || (data.org %in% c("hsa", "mmu") && net.type  %in% c("string", "innate", "huri"))){
    ppi.mat = ppi.mat[,c(2,1)];
    colnames(ppi.mat) = c("gene_id", "accession");
  }
  # fix the factor col related to library issue
  ppi.mat[,1] <- as.character(as.numeric(ppi.mat[,1]));
  ppi.mat[,2] <- as.character(ppi.mat[,2]);
  return(ppi.mat);
}

doUniprot2EntrezMapping <- function(uniprot.vec){
  
  db.map <-  queryGeneDB("entrez_uniprot", data.org);
  hit.inx <- match(uniprot.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doEntrez2UniprotMapping <- function(entrez.vec){
  # db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep="");
  # db.map <-  readRDS(db.path);
  db.map <-  queryGeneDB("entrez_uniprot", data.org);
  hit.inx <- match(entrez.vec, db.map[, "gene_id"]);
  unips <- db.map[hit.inx, "accession"];
  mode(unips) <- "character";
  return(unips);
}

doString2EntrezMapping <- function(string.vec){
  # db.path <- paste(lib.path, data.org, "/entrez_string.rds", sep="");
  # db.map <-  readRDS(db.path);
  db.map <-  queryGeneDB("entrez_string", data.org);
  hit.inx <- match(string.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doEmblGene2EntrezMapping <- function(emblgene.vec){
  # db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
  # db.map <-  readRDS(db.path);
  db.map <-  queryGeneDB("entrez_embl_gene", data.org);
  hit.inx <- match(emblgene.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}

doEmblProtein2EntrezMapping <- function(emblprotein.vec){
  # db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
  # db.map <-  readRDS(db.path);
  db.map <-  queryGeneDB("entrez_embl_protein", data.org);
  hit.inx <- match(emblprotein.vec, db.map[, "accession"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  return(entrezs);
}
