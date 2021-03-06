##################################################
## R script for NetworkAnalyst
## Description: Gene/Probe/Protein ID Annotation
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

.loadEnrLib <- function(fun.type){
  if(tolower(fun.type) == 'kegg'){ 
    current.geneset <- .loadKEGGLib();
  }else if(tolower(fun.type) == 'reactome'){ 
    current.geneset <-.loadREACTOMELib();
  }else if(tolower(fun.type) == 'motif'){ 
    current.geneset <- .loadMotifLib();
  }else{ # GO
    current.geneset <-.loadGOLib(fun.type);
  }
  return(current.geneset)
}

# Load various libaries for functional enrichment analysis
.loadKEGGLib<-function(){
  kegg.path <- paste(lib.path, data.org, "/kegg.rds", sep="");
  
  kegg.anot <- readRDS(kegg.path)
  current.setlink <- kegg.anot$link;
  current.geneset <- kegg.anot$sets;
  set.ids<- names(current.geneset); 
  names(set.ids) <- names(current.geneset) <- kegg.anot$term;
  
  current.setlink <<- current.setlink;
  current.setids <<- set.ids;
  saveRDS(current.geneset, "current_geneset.rds");
  return(current.geneset);
}

.loadREACTOMELib<-function(){
  
  reactome.path <- paste(lib.path, data.org, "/reactome.rds", sep="");
  reactome.anot <- readRDS(reactome.path)
  current.geneset <- reactome.anot$sets;
  set.ids<- names(current.geneset); 
  names(set.ids) <- names(current.geneset) <- reactome.anot$term;
  current.setlink <<- reactome.anot$link;
  current.setids <<- set.ids;
  saveRDS(current.geneset, "current_geneset.rds");
  return(current.geneset);
}

.loadMotifLib<-function(){
  
  motif.path <- paste(lib.path, data.org, "/motif_set.rds", sep="");
  motif_set<-readRDS(motif.path);
  current.geneset <- motif_set$set;
  set.ids<- names(current.geneset); 
  names(set.ids) <- names(current.geneset) <- motif_set$term;
  current.setlink <<- motif_set$link;
  current.setids <<- set.ids;
  saveRDS(current.geneset, "current_geneset.rds");
  return(current.geneset);
}

.loadGOLib<-function(onto){
  go.path <- paste(lib.path, data.org, "/go_", tolower(onto), ".rds", sep="");
  if(tolower(onto) == "bp"){
    go_bp <- readRDS(go.path);
    if(is.null(names(go_bp))){ # new go lib does not give names
      names(go_bp) <- c("link", "term", "sets");
    }
    current.link <- go_bp$link;
    current.geneset <- go_bp$sets;
    set.ids<- names(current.geneset); 
    names(set.ids) <- names(current.geneset) <- go_bp$term;
  }else if(tolower(onto) == "mf"){
    go_mf <- readRDS(go.path);
    if(is.null(names(go_mf))){
      names(go_mf) <- c("link", "term", "sets");
    }
    current.link <- go_mf$link;
    current.geneset <- go_mf$sets;
    set.ids<- names(current.geneset); 
    names(set.ids) <- names(current.geneset) <- go_mf$term;
  }else{
    go_cc <- readRDS(go.path);
    if(is.null(names(go_cc))){
      names(go_cc) <- c("link", "term", "sets");
    }
    current.link <- go_cc$link;
    current.geneset <- go_cc$sets;
    set.ids<- names(current.geneset); 
    names(set.ids) <- names(current.geneset) <- go_cc$term;
  }
  names(current.geneset) = firstup(names(current.geneset))
  names(current.geneset) = gsub("-", "_", names(current.geneset))
  names(set.ids) = firstup(names(set.ids));
  names(set.ids) = gsub("-", "_", names(set.ids))
  current.setlink <<- current.link;
  current.setids <<- set.ids;
  saveRDS(current.geneset, "current_geneset.rds");
  return(current.geneset);
}

# geneIDs is text one string, need to make to vector
performGene2ProteinMapping <- function(listNm, geneIDs, org, type){
  
  dataSet <- list();
  dataSet$orig <- geneIDs;
  current.msg <<- NULL;
  data.org <<- org;
  SetInitLib(org)
  
  listNms = vector();
  dataList <- .parseListInput(geneIDs);
  all.prot.mat <- list(); 
  for(i in 1:length(dataList)){
    dataSet$name = paste0("datalist", i);
    listNms[i] = dataSet$name;
    gene.mat <- prot.mat <- dataList[[i]];
    GeneAnotDB <-convertIdToEntrez(rownames(gene.mat), type);
    
    na.inx <- is.na(GeneAnotDB[,1]) | is.na(GeneAnotDB[,2]);
    if(sum(!na.inx) < 2){
      current.msg <<- "Less than two hits found in uniprot database. ";
      print(current.msg);
    }
    rownames(prot.mat) <- GeneAnotDB[,2];
    prot.mat <- prot.mat[!na.inx, , drop=F];
    
    # now merge duplicates
    prot.mat <- RemoveDuplicates(prot.mat, "mean", quiet=T); 
    
    if(is.null(dim(prot.mat))){
      prot.mat <- matrix(prot.mat);
    }
    
    seed.proteins <- rownames(prot.mat);
    dataSet$GeneAnotDB <- GeneAnotDB;
    dataSet$sig.mat <- gene.mat;
    dataSet$prot.mat <- prot.mat;
    dataSet$seeds.proteins <- seed.proteins;
    if(i == 1){
      all.prot.mat <- prot.mat;
      totalseed.proteins <- seed.proteins;
    }else{
      totalseed.proteins <- c(totalseed.proteins, seed.proteins);
      all.prot.mat <- rbind(all.prot.mat, prot.mat)
    }
    RegisterData(dataSet); 
  }
  all.ent.mat <<- all.prot.mat
  rownames(all.prot.mat) = doEntrez2SymbolMapping(rownames(all.prot.mat))
  all.prot.mat <<- all.prot.mat
  listNms <<- listNms
  mdata.all <- list();
  for(i in 1:length(listNms)){
    mdata.all[i] <- 1;
  }
  names(mdata.all) <- listNms;
  partialToBeSaved <<- c(partialToBeSaved, listNms);
  mdata.all <<- mdata.all
  return(totalseed.proteins)
}

doAnnotation <- function(id.vec, idType){
  if(idType %in% c("entrez", "symbol", "refseq", "genbank", "emblgene","emblprotein", "embltranscript", "orfid", "tair", "wormbase")){
    anot.id <- doGeneIDMapping(id.vec, idType);
  }else{
    anot.id <- doProbeMapping(id.vec, idType);
    names(anot.id) <- id.vec;
  }
  return(anot.id);        
}

# from probe ID to entrez ID 
doProbeMapping <- function(probe.vec, platform){
  platform.path <- paste(lib.path,  data.org, "/", platform, ".rds", sep="");
  probe.map <- readRDS(platform.path);
  if(is.null(probe.vec)){
    entrez <- probe.map[, "entrez"];
  }else{
    hit.inx <- match(probe.vec, probe.map[, "probe"]);
    entrez <- probe.map[hit.inx, "entrez"];
  }
  return(entrez);
}


queryGeneDB <- function(table.nm, data.org){
  require('RSQLite')
  conv.db <- dbConnect(SQLite(), paste(sqlite.path, data.org, "_genes.sqlite", sep="")); 
  db.map <- dbReadTable(conv.db, table.nm);
  dbDisconnect(conv.db); cleanMem();
  return(db.map)
}

# mapping between genebank, refseq and entrez
doGeneIDMapping <- function(q.vec, type){
  if(is.null(q.vec)){
    db.map <-  queryGeneDB("entrez", data.org);
    q.vec <- db.map[, "gene_id"];
    type = "entrez";
  }
  if(type == "symbol"){
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "symbol"]);
  }else if(type == "entrez"){
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
  }else{
    # note, some ID can have version number which is not in the database
    # need to strip it off NM_001402.5 => NM_001402
    q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
    q.vec <- q.mat[,1];
    if(type == "genbank"){
      db.map <-  queryGeneDB("entrez_gb", data.org);
    }else if(type == "refseq"){
      db.map <-  queryGeneDB("entrez_refseq", data.org);
    }else if(type == "emblgene"){
      db.map <-  queryGeneDB("entrez_embl_gene", data.org);
    }else if(type == "embltranscript"){
      db.map <-  queryGeneDB("entrez_embl_transcript", data.org);
    }else if(type == "emblprotein"){
      db.map <-  queryGeneDB("entrez_embl_protein", data.org);
    }else if(type == "orfid"){ # only for yeast
      db.map <-  queryGeneDB("entrez_orf", data.org);
      db.path <- paste(lib.path, data.org, "/entrez_orf.rds", sep="");
    }else if(type == "tair"){ # only for ath
      db.map <-  queryGeneDB("tair", data.org);
    }else if(type == "wormbase"){ # only for cel
      db.map <-  queryGeneDB("entrez_wormbase", data.org);
    }else{
      print("Unknown data type");
      return(0);
    }
    hit.inx <- match(q.vec, db.map[, "accession"]);
  }
  entrezs <- db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  rm(db.map, q.vec); gc();
  return(entrezs);
}

doEntrez2SymbolMapping <- function(entrez.vec){
  gene.map <-  queryGeneDB("entrez", data.org);
  gene.map[] <- lapply(gene.map, as.character)
  
  hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
  symbols <- gene.map[hit.inx, "symbol"];
  
  # if not gene symbol, use id by itself
  na.inx <- is.na(symbols);
  symbols[na.inx] <- entrez.vec[na.inx];
  return(symbols);
}

checkEntrezMatches <- function(entrez.vec){
  gene.map <-  queryGeneDB("entrez", data.org);
  gene.map[] <- lapply(gene.map, as.character)
  hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
  return(length(hit.inx));
}

doSymbol2EntrezMapping <- function(symbol.vec){
  db.map <-  queryGeneDB("entrez", data.org);
  db.map[] <- lapply(db.map, as.character);
  hit.inx <- match(symbol.vec, db.map[, "symbol"]);
  entrezs <- db.map[hit.inx, "gene_id"];
  na.inx <- is.na(entrezs);
  entrezs[na.inx] <- symbol.vec[na.inx];
  mode(entrezs) <- "character";
  return(entrezs);
}


# note, entrez.vec could contain NA/null, cannot use rownames
doEntrezIDAnot <- function(entrez.vec){
  gene.map <-  queryGeneDB("entrez", data.org);
  hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
  anot.mat <- gene.map[hit.inx, c("gene_id", "symbol", "name")];
  na.inx <- is.na(hit.inx);
  anot.mat[na.inx, "symbol"] <- entrez.vec[na.inx];
  anot.mat[na.inx, "name"] <- 'NA';
  return(anot.mat);
}

convertIdToEntrez <- function(q.vec, type){ #convert user input ids to entrez
  if(type == "entrez"){
    # need to get only our data
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
  }else if(type == "symbol"){
    db.map <-  queryGeneDB("entrez", data.org);
    hit.inx <- match(q.vec, db.map[, "symbol"]);
  }else{
    if(type == "genbank"){
      # note, some ID can have version number which is not in the database
      # need to strip it off NM_001402.5 => NM_001402
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.map <-  queryGeneDB("entrez_gb", data.org);
    }else if(type == "refseq"){
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.map <-  queryGeneDB("entrez_refseq", data.org);
    }else if(type == "emblgene"){
      db.map <-  queryGeneDB("entrez_embl_gene", data.org);
    }else if(type == "tair"){
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.map <-  queryGeneDB("tair", data.org);
    }else if(type == "embltranscript"){
      db.map <-  queryGeneDB("entrez_embl_transcript", data.org);
    }else if(type == "emblprotein"){
      db.map <-  queryGeneDB("entrez_embl_protein", data.org);
    }else if(type == "orfid"){ # only for yeast
      db.map <-  queryGeneDB("entrez_orf", data.org);
    }else if(type == "flybase"){
      db.map <-  queryGeneDB("entrez_flybase", data.org);
    }else if(type == "string"){ 
      db.map <-  queryGeneDB("entrez_string", data.org);
    }else if(type == "ecogene"){ # only for ecoli
      db.map <-  queryGeneDB("entrez_ecogene", data.org);
    }else if(type == "uniprot"){
      db.map <-  queryGeneDB("entrez_uniprot", data.org);
    }else if(type == "paelocus"){
      db.map <-  queryGeneDB("entrez", data.org);
    }else{
      print("Unknown data type");
      return(0);
    }
    hit.inx <- match(q.vec, db.map[, "accession"]);
    entrezs <- db.map[hit.inx, ];
  }
    entrezs <- db.map[hit.inx, ]; 
  if(type == "entrez"){
    entrezs = entrezs[,c(1,1)];
  }else{
    entrezs = entrezs[,c(2,1)];
  }
  na.inx <- is.na(entrezs[,1]);
  entrezs[,1][na.inx] <- q.vec[na.inx];
  na.inx <- is.na(entrezs[,2]);
  entrezs[,2][na.inx] <- q.vec[na.inx];
  colnames(entrezs) <- c("accession", "gene_id")
  return(entrezs);
}

PerformListAnnot <- function(listNm, org, geneIDs, type){
  dataSet <- list();
  dataSet$orig <- "";
  current.msg <<- NULL;
  data.org <<- org;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  SetInitLib(org)
  listNms <- multiFileNamesU;
  numOfLists <<-length(multiFileNamesU);
  notOk = 0
  for(i in 1:length(listNms)){
    dataSet = readRDS(listNms[i])
    dataSet$name = listNms[i];
    gene.mat <- prot.mat <- dataSet$prot.mat;
    GeneAnotDB <-convertIdToEntrez(rownames(gene.mat), type);
    
    na.inx <- is.na(GeneAnotDB[,1]) | is.na(GeneAnotDB[,2]);
    if(sum(!na.inx) < 2){
      current.msg <<- paste0("Less than two hits found in database for ", listNms[i]);
      print(current.msg);
      return(0);
    }
    rownames(prot.mat) <- GeneAnotDB[,2];
    prot.mat <- prot.mat[!na.inx, , drop=F];
    
    # now merge duplicates
    listInxU <<- listNms[i];
    prot.mat <- RemoveDuplicates(prot.mat, "mean", quiet=F); 
    
    if(is.null(dim(prot.mat))){
      prot.mat <- matrix(prot.mat);
    }
    seed.proteins <- rownames(prot.mat);
    dataSet$GeneAnotDB <- GeneAnotDB;
    dataSet$sig.mat <- gene.mat;
    dataSet$prot.mat <- prot.mat;
    dataSet$seeds.proteins <- seed.proteins;
    if(i == 1){
      all.prot.mat <- prot.mat;
      totalseed.proteins = seed.proteins
    }else{
      totalseed.proteins  = c(totalseed.proteins, seed.proteins);
      all.prot.mat <- rbind(all.prot.mat, prot.mat)
    }
    RegisterData(dataSet); 
  }
  all.ent.mat <<- all.prot.mat;
  rownames(all.prot.mat) = doEntrez2SymbolMapping(rownames(all.prot.mat))
  all.prot.mat <<- all.prot.mat;
  listNms <<- listNms
  mdata.all <- list();
  for(i in 1:length(listNms)){
    mdata.all[i] <- 1;
  }
  names(mdata.all) <- listNms;
  mdata.all <<- mdata.all
  return(totalseed.proteins)
}
