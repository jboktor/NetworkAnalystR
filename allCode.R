##################################################
## R script for NetworkAnalyst
## Description: Gene/Probe/Protein ID Annotation
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

LoadEnrLib <- function(fun.type){
  if(tolower(fun.type) == 'kegg'){ 
    current.geneset <- LoadKEGGLib();
  }else if(tolower(fun.type) == 'reactome'){ 
    current.geneset <-LoadREACTOMELib();
  }else if(tolower(fun.type) == 'motif'){ 
    current.geneset <- LoadMotifLib();
  }else{ # GO
    current.geneset <-LoadGOLib(fun.type);
  }
  return(current.geneset)
}

# Load various libaries for functional enrichment analysis
LoadKEGGLib<-function(){
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

LoadREACTOMELib<-function(){
  
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

LoadMotifLib<-function(){
  
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

LoadGOLib<-function(onto){
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
      totalseed.proteins = seed.proteins
    }else{
      totalseed.proteins  = c(totalseed.proteins, seed.proteins);
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
  mdata.all <<- mdata.all
  return(totalseed.proteins)
}

GetNumOfLists <- function(){
  return(numOfLists)
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
  
  conv.db <- dbConnect(SQLite(), paste(genesdb.path, data.org, "_genes.sqlite", sep="")); 
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
  entrezs=db.map[hit.inx, "gene_id"];
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
  colnames(entrezs) <- c("accession", "gene_id")
  return(entrezs);
}

# note: hit.query, resTable must synchronize
# ora.vec should contains entrez ids, named by their gene symbols
PerformEnrichAnalysis <- function(file.nm, fun.type, ora.vec){
  # prepare lib
  if(tolower(fun.type) == 'kegg'){ 
    current.geneset <-LoadKEGGLib();
  }else if(tolower(fun.type) == 'reactome'){ 
    current.geneset <-LoadREACTOMELib();
  }else if(tolower(fun.type) == 'motif'){ 
    current.geneset <-LoadMotifLib();
  }else{ # GO
    current.geneset <- LoadGOLib(fun.type);
  }
  
  # prepare query
  ora.nms <- names(ora.vec);
  
  # prepare for the result table
  set.size<-length(current.geneset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.geneset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "P.Value", "FDR");
  
  # need to cut to the universe covered by the pathways, not all genes
  current.universe <- unique(unlist(current.geneset)); 
  hits.inx <- ora.vec %in% current.universe;
  ora.vec <- ora.vec[hits.inx];
  ora.nms <- ora.nms[hits.inx];
  
  q.size<-length(ora.vec);
  
  # get the matched query for each pathway
  hits.query <- lapply(current.geneset, 
                       function(x) {
                         ora.nms[ora.vec%in%unlist(x)];
                       }
  );
  
  saveRDS(hits.query, "hits_query.rds");
  
  names(hits.query) <- names(current.geneset);
  hit.num<-unlist(lapply(hits.query, function(x){length(unique(x))}), use.names=FALSE);
  
  # total unique gene number
  uniq.count <- length(current.universe);
  
  # unique gene count in each pathway
  set.size <- unlist(lapply(current.geneset, length));
  
  res.mat[,1]<-set.size;
  res.mat[,2]<-q.size*(set.size/uniq.count);
  res.mat[,3]<-hit.num;
  
  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
  res.mat[,4]<- raw.pvals;
  res.mat[,5] <- p.adjust(raw.pvals, "fdr");
  
  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];
  
  if(nrow(res.mat)> 1){
    # order by p value
    ord.inx<-order(res.mat[,4]);
    res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];
    
    imp.inx <- res.mat[,4] <= 0.05;
    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      if(sum(imp.inx) > 120){
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }else{
    return(0);
  }
  
  #get gene symbols
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
  res.mat[,"Hits"] = res.mat[,"Hits"]
  enr.mat <<- res.mat
  current.msg <<- "Functional enrichment analysis was completed";
  
  # write json
  require(RJSONIO);
  fun.anot = hits.query; 
  total = resTable[,2]; if(length(total) ==1) { total <- matrix(total) };
  fun.pval = resTable[,5]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  fun.padj = resTable[,6]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  hit.num = resTable[,4]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
  fun.ids <- as.vector(current.setids[names(fun.anot)]); 
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  json.res <- list(
    fun.link = current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    hit.num = hit.num,
    total= total
  );
  json.mat <- toJSON(json.res, .na='null');
  json.nm <- paste(file.nm, ".json", sep="");
  
  sink(json.nm)
  cat(json.mat);
  sink();
  
  # write csv
  fun.hits <<- hits.query;
  fun.pval <<- resTable[,5];
  hit.num <<- resTable[,4];
  csv.nm <- paste(file.nm, ".csv", sep="");    
  write.csv(resTable, file=csv.nm, row.names=F);
  
  return(1);
}

SetMultiNames <-function(names){
  multiNamesU <<- names;
}

ReadListData <- function(fileName) {
  dat1 <- data.table::fread(fileName, header=FALSE, check.names=FALSE, data.table=FALSE);
  dataSet$name = fileName
  rowNms = dat1[,1]
  if(length(dat1) == 1){
    dat1[,1] = 0
  }else{
    dat1[,1] = dat1[,2]
    dat1 = dat1[,-2];
  }
  dataSet$prot.mat = as.matrix(dat1)
  rownames(dataSet$prot.mat) = rowNms;
  saveRDS(dataSet, file=fileName); # keep original copy, not in mem
  return(1)
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
  all.ent.mat <<- all.prot.mat
  rownames(all.prot.mat) = doEntrez2SymbolMapping(rownames(all.prot.mat))
  all.prot.mat <<- all.prot.mat
  listNms <<- listNms
  mdata.all <- list();
  for(i in 1:length(listNms)){
    mdata.all[i] <- 1;
  }
  names(mdata.all) <- listNms;
  mdata.all <<- mdata.all
  return(totalseed.proteins)
}

##################################################
## R script for NetworkAnalyst
## Description: prepare data for chord diagram
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

PrepareChordDataInit<-function(){
  # create a list store all possible combination (for a max of 4)
  # note, the whole region are divided into small areas (16 for 4; 7 for 3, 3 for 2)
  # for instance:
  # a: a unique (no b, no c)
  # ab: a and b, no c
  newDat <- list();
  if(anal.type == "metadata"){
    sel.inx <- mdata.all==1;
    sel.nms <- names(mdata.all)[sel.inx];
  }else{
    sel.nms <- listNms;
  }
  venn.genenb <- vector()
  sel.dats <- list();
  for(i in 1:length(sel.nms)){
    nm = sel.nms[i]
    dataSet <- readRDS(nm);
    if(anal.type == "metadata"){
      sel.dats[[nm]] <- rownames(dataSet$sig.mat)
    }else{
      sel.dats[[nm]] <- rownames(dataSet$prot.mat)
    }
    venn.genenb[i] = length(sel.dats[[nm]])
  }
  if(anal.type == "metadata" & meta.selected){
    sel.dats[["meta_dat"]] <- as.character(meta.stat$de);
    venn.genenb[length(venn.genenb) + 1] = length(as.character(meta.stat$de))
  }
  
  chord.list <<- sel.dats;
  chord.genenb <<- venn.genenb
  chord.list.up <<- sel.dats;
  chord.genenb.up <<- venn.genenb
  return(PrepareChordData());
}

PrepareSelChordData<-function(selectedNms){
  newDat <- list();
  sel.nms <- unlist(strsplit(selectedNms, ";"));
  nm.vec <<- sel.nms;
  SelectData();
  venn.genenb <- vector()
  sel.dats <- list();
  for(i in 1:length(sel.nms)){
    nm = sel.nms[i]
    if(nm != "meta_dat"){
      dataSet <- readRDS(nm);
      if(anal.type == "metadata"){
        sel.dats[[nm]] <- rownames(dataSet$sig.mat)
      }else{
        sel.dats[[nm]] <- rownames(dataSet$prot.mat)
      }
      venn.genenb[i] = length(sel.dats[[nm]])
    }else{
      sel.dats[[nm]] <- as.character(meta.stat$de);
      venn.genenb[i] = length(as.character(meta.stat$de))
    }
  }
  
  chord.list.up <<- sel.dats;
  chord.genenb.up <<- venn.genenb
  return(PrepareChordData());
}

GetSelectedDataNamesChord <- function(){
  return(paste(names(chord.list), collapse=";"));
}

GetSelectedDataGeneNumberChord<- function(){
  return(paste(chord.genenb, collapse=";"));
}

GetSelectedDataNamesUpdatedChord <- function(){
  return(paste(names(chord.list.up), collapse=";"));
}

GetSelectedDataGeneNumberUpdatedChord <- function(){
  return(paste(chord.genenb.up, collapse=";"));
}

GetChordFileCount <- function(){
  return(chord_count-1);
}

PrepareChordData <-function(){
  if(anal.type == "metadata"){
    res <- PrepareMetaChordData();
  }else{
    res <- PrepareListChordData();
  }
  
  if(is.null(res)){
    return(0);
  }
  require(RJSONIO);
  chordData <- res$chordData;
  fileNm = paste0("networkanalyst_chorddata_",chord_count,".json")
  sink(fileNm);
  cat(toJSON(chordData));
  sink();
  
  lookup <- res$lookup;
  fileNm2 = paste0("networkanalyst_chord_lookup_",chord_count,".json")
  sink(fileNm2);
  cat(toJSON(lookup));
  chord_count <<- chord_count + 1
  sink();
  return(1);
}

PrepareChordDataFromList <- function(newDat, uniq.enIDs){
  # now combine lists into a single matrix
  all.nms <- listNms
  nm.used <- NULL;
  for(i in 1:length(all.nms)){
    dataNm <- all.nms[i];
    if(mdata.all[[dataNm]] == 1){ # selected
      nm.used <- c(nm.used, dataNm);
    }
  }
  if(anal.type == "metadata" & meta.selected){
    nm.used <- c(nm.used, "meta");
  }
  hit.mat <- fc.mat <- matrix(0, nrow=length(nm.used), ncol=length(uniq.enIDs));
  nm.mat <- matrix(NA, nrow=length(nm.used), ncol=length(uniq.enIDs));
  shared.ids <- uniq.enIDs;
  shared.sbls <- doEntrez2SymbolMapping(uniq.enIDs);
  
  dataNms <-nm.used;
  colnames(hit.mat) <- colnames(fc.mat) <- colnames(nm.mat) <- shared.sbls;
  rownames(hit.mat) <- rownames(fc.mat) <- rownames(nm.mat) <- dataNms;
  names(shared.sbls) <- shared.ids;
  
  # push to parent env.
  shared.sbls <<- shared.sbls;
  
  for(i in 1:length(dataNms)){
    nm <- dataNms[i];
    exp.vec <- newDat[[nm]];
    hit.inx <- shared.ids %in% names(exp.vec);
    hit.mat[i, hit.inx] <- 1;
    nm.mat[i, hit.inx] <- nm;
    fc.mat[i, hit.inx] <- as.numeric(exp.vec);
  }
  
  ##############Links#######
  ## links & arcs data 
  ###########################
  
  links <- arcs <- list();
  
  # arcs data
  de.num <- apply(hit.mat, 1, sum);
  de.prct <- as.numeric(de.num/sum(de.num));
  
  for(i in 1:length(dataNms)){
    arcs[[i]] <- list(
      name = dataNms[i],
      label = dataNms[i],
      size = as.numeric(de.num[i]),
      val = de.prct[i]
    )
  }
  
  # now, get unique names by each inserting in the gene names
  lnk.genes <- colnames(nm.mat);
  nm.mat <- apply(nm.mat, 1, 
                  function(x){
                    gd.hits <- !is.na(x); 
                    x[gd.hits] <-paste(x[gd.hits], "*", lnk.genes[gd.hits], sep="");
                    x;
                  });
  
  nm.mat <- unname(t(nm.mat));
  for(l in 1:nrow(nm.mat)){
    x <- nm.mat[l,];      
    hit.inx <- !is.na(x);
    links$name <- c(links$name, x[hit.inx]);
    
    y.mat <- nm.mat[,hit.inx,drop=F];
    new.mat <- apply(y.mat, 2, 
                     function(d){
                       myd <- d[-l];
                       if(all(is.na(myd))){ # link to itself?
                         d[l];
                       }else{
                         myd[!is.na(myd)];
                       }
                     });
    links$imports <- c(links$imports,new.mat);
  }
  
  # need to reformat
  links.new <- list();
  for(l in 1:length(links$name)){
    impt<-links$imports[[l]];
    if(length(impt) == 1){
      impt <- as.matrix(impt);
    }
    links.new[[l]] <- list(
      name=links$name[l],
      imports=impt
    )
  }
  
  ###################
  #colors & labels
  require(RColorBrewer);  
  if(length(dataNms) > 9){
    colors <- brewer.pal(length(dataNms),"Set3");
  }else if (length(dataNms) > 2){
    colors <- brewer.pal(length(dataNms),"Dark2");
  }else {
    colors <- c("#7570B3", "#D95F02");
  }
  labels <- rep("", length(dataNms));
  names(labels) <- names(colors) <- dataNms;
  
  ###################
  # links done!, now get weights/pvalues/data for every pairs of chord
  dataNms <- rownames(fc.mat);
  geneNms <- colnames(fc.mat);
  entrezIDs <- names(shared.sbls);
  weights <- data <- list();
  
  # also need to setup lookup json
  lookup <- list();
  for(i in 1:length(dataNms)){
    nlst <- vector(length = length(dataNms)-1, mode = "list"); 
    nm <- dataNms[i];
    names(nlst) <- dataNms[-i]
    lookup[[nm]] <- nlst;
  }
  
  item.count <- 0;
  for(i in 1:ncol(nm.mat)){
    geneNm <- geneNms[i]; 
    entrez <- entrezIDs[i];
    hit.inx <- !is.na(nm.mat[,i]);
    fc.vals <- as.numeric(fc.mat[hit.inx,i]);
    
    datNms <- dataNms[hit.inx];
    
    if(length(datNms) == 1){
      item.count <- item.count + 1;
      id <-paste(datNms,"*",datNms,"*",geneNm, sep="")
      #weights[[id]]<-fc.vals;
      weights[[id]]<-1;
      lookup[[datNms]][[datNms]][[geneNm]] <- list(
        entrez = entrez,
        genesymbol = geneNm,
        weight = fc.vals
      )
      data[[item.count]]<-list(
        fromModule = datNms,
        toModule= datNms,
        arc=geneNm,
        weight=fc.vals,
        fromColor=as.character(colors[datNms]),
        toColor=as.character(colors[datNms])
      );
    }else{
      for(m in 1:(length(datNms)-1)){
        nm1 <- datNms[m];
        fc1 <- fc.vals[m];
        for(n in (m+1):length(datNms)){
          item.count <- item.count + 1;
          nm2 <- datNms[n];
          fc2 <- fc.vals[n];
          id1<- paste(nm1,"*",nm2,"*",geneNm, sep="");
          id2<- paste(nm2,"*",nm1,"*",geneNm, sep="");
          wt <- 1;
          weights[[id1]]<-weights[[id2]]<-1;
          
          lookup[[nm1]][[nm2]][[geneNm]] <- list(
            entrez = entrez,
            genesymbol = geneNm,
            weight = 1
          )
          lookup[[nm2]][[nm1]][[geneNm]] <- list(
            entrez = entrez,
            genesymbol = geneNm,
            weight = 1
          )
          # only from large module (src) to small (no need for reverse)
          if(de.num[nm1] >= de.num[nm2]){
            data[[item.count]]<-list(
              fromModule = nm1,
              toModule= nm2,
              arc=geneNm,
              weight=wt,
              fromColor=as.character(colors[nm1]),
              toColor=as.character(colors[nm2])
            );
          }else{
            data[[item.count]]<-list(
              fromModule = nm2,
              toModule=nm1,
              arc=geneNm,
              weight=wt,
              fromColor=as.character(colors[nm2]),
              toColor=as.character(colors[nm1])
            );
          }
        }
      }
    }
  }
  
  # need to remove 3rd level list name 
  for(m in 1:length(dataNms)){ # need to remove itself
    nm1 <- dataNms[m];
    for(n in m:length(dataNms)){
      nm2 <- dataNms[n];
      nlist <- lookup[[nm1]][[nm2]];
      lookup[[nm1]][[nm2]] <- lookup[[nm2]][[nm1]] <- unname(nlist);
    }
  }
  
  chordData <- list(
    arcs = arcs,
    links = links.new,
    weights = weights,
    data = data,
    colors = colors,
    labels = labels
  )
  
  return (
    list(
      chordData = chordData,
      lookup = lookup
    )
  );
}

PrepareMetaChordData <- function(){
  BHth <- GlobalCutOff$BHth;
  logFC <- GlobalCutOff$logFC;
  sel.nms <- names(mdata.all)[mdata.all==1];
  row.nm <- length(sel.nms);
  
  # update inmex.ind for only selected dataset
  hit.inx <-  names(inmex.ind) %in% sel.nms;
  inmex.ind <- inmex.ind[hit.inx];
  dataNms <- names(inmex.ind);
  newNms <- substring(dataNms,0, nchar(dataNms)-4);
  
  if(meta.selected){
    newNms <- c(newNms, "meta");
    row.nm <- row.nm + 1;
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  hit.mat <- matrix(0, nrow=row.nm, ncol=nrow(inmex.meta$data));
  nm.mat <- matrix(NA, nrow=row.nm, ncol=nrow(inmex.meta$data));
  shared.ids <- rownames(inmex.meta$data);
  shared.sbls <- inmex.meta$gene.symbls;
  
  colnames(hit.mat) <- colnames(nm.mat) <- shared.sbls;
  rownames(hit.mat) <- rownames(nm.mat) <- newNms;
  
  chord.vec.nms <- newNms;
  chord.res <- list();
  tot.count <- 0; # need to see how many chords
  for(i in 1:length(inmex.ind)){
    dataSet = readRDS(sel.nms[i]);
    nm <- newNms[i];
    res.mat <- inmex.ind[[i]][shared.ids, ];
    hit.inx <- abs(res.mat[,1]) >= logFC & res.mat[,2] <= dataSet$pval;
    hit.mat[i, hit.inx] <- 1;
    nm.mat[i, hit.inx] <- nm;
    chord.res[[nm]] <- shared.sbls[hit.inx];
    tot.count <- tot.count + sum(hit.inx);
  }
  
  if(meta.selected){# add meta
    i  <- length(inmex.ind) + 1;
    hit.inx <- shared.ids %in% meta.stat$de;
    hit.mat[i, hit.inx] <- 1;
    nm.mat[i, hit.inx] <- newNms[i];
    chord.res[[newNms[i]]] <- shared.sbls[hit.inx];
    tot.count <- tot.count + sum(hit.inx);
  }
  
  if(tot.count > 2000){
    current.msg <<- paste("Chord diagrams is effective to display relationships for less than 1000 genes. The results contain", tot.count, 
                          "of genes (max. allowed: 2000). You can try Venn diagram or adjust threshold to select most significant genes.")
    print(current.msg);
    return(NULL);
  }
  
  # keep gene with at least 1 hit across all datasets
  keep.inx <- apply(hit.mat, 2, sum) > 0;
  hit.mat <- hit.mat[,keep.inx];
  nm.mat <- nm.mat[,keep.inx];
  shared.sbls <- shared.sbls[keep.inx];
  
  ##############Links#######
  #links & arcs data 
  links <- arcs <- list();
  
  # arcs data
  de.num <- apply(hit.mat, 1, sum);
  de.prct <- as.numeric(de.num/sum(de.num));
  
  for(i in 1:length(newNms)){
    arcs[[i]] <- list(
      name = newNms[i],
      label = newNms[i],
      size = as.numeric(de.num[i]),
      val = de.prct[i]
    )
  }
  
  # now, get unique names by each inserting in the gene names
  lnk.genes <- colnames(nm.mat);
  nm.mat <- apply(nm.mat, 1, 
                  function(x){
                    gd.hits <- !is.na(x); 
                    x[gd.hits] <-paste(x[gd.hits], "*", lnk.genes[gd.hits], sep="");
                    x;
                  });
  
  nm.mat <- unname(t(nm.mat));
  for(l in 1:nrow(nm.mat)){
    x <- nm.mat[l,];      
    hit.inx <- !is.na(x);
    links$name <- c(links$name, x[hit.inx]);
    
    y.mat <- nm.mat[,hit.inx,drop=F];
    new.mat <- apply(y.mat, 2, 
                     function(d){
                       myd <- d[-l];
                       if(all(is.na(myd))){ # link to itself?
                         d[l];
                       }else{
                         myd[!is.na(myd)];
                       }
                     });
    links$imports <- c(links$imports,new.mat);
  }
  
  # need to reformat
  links.new <- list();
  for(l in 1:length(links$name)){
    impt<-links$imports[[l]];
    if(length(impt) == 1){
      impt <- as.matrix(impt);
    }
    links.new[[l]] <- list(
      name=links$name[l],
      imports=impt
    )
  }
  
  ###################
  #colors & labels
  require(RColorBrewer);  
  if(length(newNms) > 9){
    colors <- brewer.pal(length(newNms),"Set3");
  }else if (length(newNms) > 2){
    colors <- brewer.pal(length(newNms),"Dark2");
  }else {
    colors <- c("#7570B3", "#D95F02");
  }
  labels <- rep("", length(newNms));
  names(labels) <- names(colors) <- newNms;
  
  ###################
  # links done!, now get pvalues/data for every pairs of chord
  dataNms <- rownames(hit.mat);
  geneNms <- colnames(hit.mat);
  entrezIDs <- names(shared.sbls);
  weights <- data <- list();
  
  # also need to setup lookup json
  lookup <- list();
  for(i in 1:length(dataNms)){
    nlst <- vector(length = length(dataNms)-1, mode = "list"); 
    nm <- dataNms[i];
    names(nlst) <- dataNms[-i]
    lookup[[nm]] <- nlst;
  }
  
  item.count <- 0;
  for(i in 1:ncol(nm.mat)){
    geneNm <- geneNms[i]; 
    entrez <- entrezIDs[i];
    hit.inx <- !is.na(nm.mat[,i]);
    datNms <- dataNms[hit.inx];
    
    if(length(datNms) == 1){
      item.count <- item.count + 1;
      id <-paste(datNms,"*",datNms,"*",geneNm, sep="");
      weights[[id]]<-1;
      lookup[[datNms]][[datNms]][[geneNm]] <- list(
        entrez = entrez,
        genesymbol = geneNm
      )
      data[[item.count]]<-list(
        fromModule = datNms,
        toModule= datNms,
        arc=geneNm,
        fromColor=as.character(colors[datNms]),
        toColor=as.character(colors[datNms])
      );
    }else{
      for(m in 1:(length(datNms)-1)){
        nm1 <- datNms[m];
        for(n in (m+1):length(datNms)){
          item.count <- item.count + 1;
          nm2 <- datNms[n];
          id1<- paste(nm1,"*",nm2,"*",geneNm, sep="");
          id2<- paste(nm2,"*",nm1,"*",geneNm, sep="");
          weights[[id1]]<-weights[[id2]]<-1;
          
          lookup[[nm1]][[nm2]][[geneNm]] <- list(
            entrez = entrez,
            genesymbol = geneNm
          )
          lookup[[nm2]][[nm1]][[geneNm]] <- list(
            entrez = entrez,
            genesymbol = geneNm
          )
          # only from large module (src) to small (no need for reverse)
          if(de.num[nm1] >= de.num[nm2]){
            data[[item.count]]<-list(
              fromModule = nm1,
              toModule= nm2,
              arc=geneNm,
              fromColor=as.character(colors[nm1]),
              toColor=as.character(colors[nm2])
            );
          }else{
            data[[item.count]]<-list(
              fromModule = nm2,
              toModule=nm1,
              arc=geneNm,
              fromColor=as.character(colors[nm2]),
              toColor=as.character(colors[nm1])
            );
          }
        }
      }
    }
  }
  
  # need to remove 3rd level list name 
  for(m in 1:length(dataNms)){ # need to remove itself
    nm1 <- dataNms[m];
    for(n in m:length(dataNms)){
      nm2 <- dataNms[n];
      nlist <- lookup[[nm1]][[nm2]];
      lookup[[nm1]][[nm2]] <- lookup[[nm2]][[nm1]] <- unname(nlist);
    }
  }
  
  chordData <- list(
    arcs = arcs,
    links = links.new,
    data = data,
    weights = weights,
    colors = colors,
    labels = labels
  )
  
  chord.vec.nms <<- chord.vec.nms;
  chord.res <<- chord.res;
  
  return (
    list(
      chordData = chordData,
      lookup = lookup
    )
  );
}

PerformChordEnrichment <- function(file.nm, fun.type, IDs){
  gene.vec <- unlist(strsplit(IDs, "; "));
  sym.vec <- doEntrez2SymbolMapping(gene.vec);
  names(gene.vec) <- sym.vec;
  res <- PerformEnrichAnalysis(file.nm, fun.type, gene.vec);
  return(res);
}

PrepareListChordData <- function(){
  
  all.enIDs <- NULL;
  newDat <- list();
  tot.count <- 0;
  all.nms <- names(mdata.all)[mdata.all==1];
  for(i in 1:length(all.nms)){
    dataNm <- all.nms[i];
    dataSet <- readRDS(dataNm);
    gene.mat <- dataSet$prot.mat;
    
    # convert to entrez
    expr.val <- gene.mat[,1];
    en.ids <- rownames(gene.mat);
    
    names(expr.val) <- en.ids;
    newDat[[dataNm]] <- expr.val;
    
    all.enIDs <- c(all.enIDs, en.ids);
    tot.count <- tot.count + nrow(gene.mat);
    
    if(tot.count > 2000){
      current.msg <<- paste("Chord diagrams is effective to display relationships for less than 1000 items. The results contain", tot.count, 
                            "of genes (max. allowed: 2000). You can try Venn diagram instead.")
      return(NULL);
    }
  }
  PrepareChordDataFromList(newDat, unique(all.enIDs));
}

CalculateDEgeneSet <- function(nms, operation, refNm, filenm){
  nms <- strsplit(nms, ";")[[1]];
  if(anal.type == "metadata"){
    com.smbls <- PerformSetOperation_Data(nms, operation, refNm);
  }else{
    com.smbls <- PerformSetOperation_List(nms, operation, refNm);
  }
  
  sink(filenm);
  cat(toJSON(com.smbls));
  sink();
}

PerformSetOperation_List <- function(nms, operation, refNm){
  all.nms <- names(mdata.all);
  include.inx <- all.nms %in% nms;
  my.vec <- all.nms[include.inx];
  if(operation == "diff"){ # make sure reference is the first
    inx <- which(my.vec == refNm);
    my.vec <- my.vec[-inx];
  }
  com.ids <- NULL;
  ids.list <- list()
  for(i in 1:length(my.vec)){
    dataSet <- readRDS(my.vec[i]);
    if(operation == "diff"){
      ids.list[[i]]=dataSet$GeneAnotDB[,"gene_id"]
      #com.ids <- setdiff(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
    }else if(is.null(com.ids)){
      com.ids <- dataSet$GeneAnotDB[,"gene_id"];
    }else{
      if(operation == "intersect"){
        com.ids <- intersect(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
      }else if(operation == "union"){
        com.ids <- union(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
      }
    }
  }
  if(operation == "diff"){
    dataSet <- readRDS(refNm);
    ids <- unique(unlist(ids.list));
    com.ids <-setdiff(dataSet$GeneAnotDB[,"gene_id"], ids);
  }
  
  com.ids <- as.character(com.ids[!is.na(com.ids)]); # make sure it is interpreted as name not index
  com.symbols <- doEntrez2SymbolMapping(com.ids);
  names(com.symbols) = com.ids
  
  com.symbols<-com.symbols[!is.null(com.symbols)];
  return(com.symbols);
}

PerformSetOperation_Data <- function(nms, operation, refNm){
  include.inx <- chord.vec.nms %in% nms;
  my.vec <- chord.vec.nms[include.inx];
  my.vec <- paste0(my.vec,".txt")
  refNm = paste0(refNm, ".txt")
  if(operation == "diff"){ # make sure reference is the first
    inx <- which(my.vec == refNm);
    my.vec <- my.vec[-inx];
  }
  com.ids <- NULL;
  ids.list <- list()
  for(nm in my.vec){
    dataSet = readRDS(nm);
    if(operation == "diff"){
      ids.list[[nm]]=rownames(dataSet$sig.mat);
      #com.ids <- setdiff(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
    }else if(is.null(com.ids)){
      com.ids <- rownames(dataSet$sig.mat);
    }else{
      if(operation == "intersect"){
        com.ids <- intersect(com.ids, rownames(dataSet$sig.mat));
      }else if(operation=="union"){
        com.ids <- union(com.ids, rownames(dataSet$sig.mat));
      }
    }
  }
  if(operation == "diff"){
    dataSet <- readRDS(refNm);
    ids <- unique(unlist(ids.list));
    com.ids <-setdiff(rownames(dataSet$sig.mat), ids);
  } 
  com.symbols <- doEntrez2SymbolMapping(com.ids);
  names(com.symbols) = com.ids;
  return(com.symbols);
}

##################################################
## R script for NetworkAnalyst
## Description: functions only for single gene expression data
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# read tab delimited file
# stored in dataSet list object
# can have many classes, stored in meta.info  
ReadTabExpressData <- function(fileName) {
  
  dataSet <- .readTabData(fileName);
  
  # rename data to data.orig
  int.mat <- dataSet$data;
  dataSet$data <- NULL;
  
  msg <- paste("a total of ", ncol(int.mat), " samples and ", nrow(int.mat), " features were found");
  
  # remove NA, null
  row.nas <- apply(is.na(int.mat)|is.null(int.mat), 1, sum);
  good.inx<- row.nas/ncol(int.mat) < 0.5;
  if(sum(!good.inx) > 0){
    int.mat <- int.mat[good.inx,];
    msg <- c(msg, paste("removed ", sum(!good.inx), " features with over 50% missing values"));
  }
  # remove constant values
  filter.val <- apply(int.mat, 1, IQR, na.rm=T);
  good.inx2 <- filter.val > 0;
  if(sum(!good.inx2) > 0){
    int.mat <- int.mat[good.inx2,];
    msg <- c(msg, paste("removed ", sum(!good.inx2), " features with constant values"));
  }
  
  if(nrow(int.mat) > 5000){
    filter.val <- filter.val[good.inx2];
    rk <- rank(-filter.val, ties.method='random');
    
    var.num <- nrow(int.mat);
    kept.num <- 0.95*var.num;
    int.mat <- int.mat[rk < kept.num, ];
    # msg <- c(msg, paste("removed 5% features with near-constant values"));
  }
  
  minVal <- min(int.mat, na.rm=T);
  na.inx <- is.na(int.mat);
  if(sum(na.inx) > 0){
    int.mat[na.inx] <- minVal/2;
    # msg <- c(msg, "the remaining", sum(na.inx), "missing variables were replaced with data min");
  }
  current.msg <<- paste(msg, collapse="; ");
  data.proc <- RemoveDuplicates(int.mat, "mean", quiet=T);
  dataSet$smpl.num <- ncol(data.proc);
  
  # save processed data for download user option
  write.csv(data.proc, file="data_processed.csv");
  saveRDS(data.proc, "data.proc.rds");
  
  dataSet <<- dataSet;
  return (1);
}

GetSampleNumber <-function(){
  return(dataSet$smpl.num);
}

GetMetaInfo <- function(){
  return(colnames(dataSet$meta.info));
}

# read in the data and perform
# gene ID mapping using built in libraries
# matchMin is minimal matched probe (%)
# idType: INVEX supported ID types
# lvlOpt: "NA" to keep original, other values will merge original ID to entrez gene IDs

# return the total matched gene number
# note: unmapped IDs will be retained as 
# original label (i.e. intergenic regions) in further analysis

PerformDataAnnot <- function(org, dataType, idType, lvlOpt){
  data.org <<- org;
  SetInitLib(org)
  
  dataSet$type <- dataType;
  dataSet$id.orig <- dataSet$id.current <- idType;
  dataSet$annotated <- F;
  # should not contain duplicates, however sanity check
  data.proc <- readRDS("data.proc.rds");
  dataSet$data.anot <- data.proc;
  
  if (org != 'NA' & idType != 'NA'){
    feature.vec <- rownames(data.proc);
    anot.id <- doAnnotation(feature.vec, idType);
    
    #dataSet$annotation <- anot.id; 
    saveRDS(anot.id, "annotation.rds");
    
    hit.inx <- !is.na(anot.id);
    matched.len <- sum(hit.inx);
    perct <- round(matched.len/length(feature.vec),3)*100;
    thresh <- 0.1 # previous value of 0.25 is causing challenges 
    #for datasets like Ppromelas with low annotation quality
    if (matched.len < length(feature.vec)*thresh){
      current.msg <<- paste('Only ', perct, '% ID were matched. You may want to choose another ID type or use default.', sep=""); 
    } else {
      current.msg <<- paste("ID annotation: ", "Total [", length(anot.id), 
                            "] Matched [", matched.len, "] Unmatched [", sum(!hit.inx),"]", collapse="\n");    
      
      if (lvlOpt != 'NA' | idType == "entrez"){
        # do actual summarization to gene level
        
        matched.entrez <- anot.id[hit.inx];
        data.anot <- data.proc[hit.inx,];
        rownames(data.anot) <- matched.entrez;
        current.msg <<- paste(current.msg, "Data is now transformed to gene-level (Entrez) expression.");
        
        dataSet$data.anot <- RemoveDuplicates(data.anot, lvlOpt, quiet=F);
        dataSet$id.current <- "entrez";
        dataSet$annotated <- T; 
      } else {
        current.msg <<- paste(current.msg, "No gene level summarization was performed.");
      }
    }
  } else { # no conversion will be performed
    feature.vec <- rownames(data.proc);
    anot.id = feature.vec
    perct <- 100;
    hit.inx <- !is.na(anot.id);
    matched.len <- length(feature.vec); # dummies
    minLvl <- 1;
    current.msg <<- paste("No annotation was performed. Make sure organism and gene ID are specified correctly!"); 
  }
  # need to save the ids (mixed gene annotation and original id) 
  # in case, users needs to keep unannotated features
  # this need to be updated to gether with data from now on
  dataSet$data.norm <- dataSet$data.anot;
  dataSet <<- dataSet;
  
  saveRDS(dataSet$data.anot, file="orig.data.anot"); # keep original copy, not in mem
  
  totalCount =  sum(colSums(dataSet$data.anot));
  avgCount = sum(colSums(dataSet$data.anot))/ ncol(dataSet$data.anot);
  minCount = min(colSums(dataSet$data.anot))
  maxCount = max(colSums(dataSet$data.anot))
  
  if(length(dataSet$meta.info)==1){
    lvls = paste(levels(dataSet$meta.info[,1]),collapse="; ")
  }else{
    conc1 = paste0("<b>", colnames(dataSet$meta.info)[1], "</b>", ": ", paste(levels(dataSet$meta.info[,1]), collapse="; "))
    conc2 = paste0("<b>", colnames(dataSet$meta.info)[2], "</b>", ": ", paste(levels(dataSet$meta.info[,2]), collapse="; "))
    lvls = paste("Two factors found -", conc1, conc2)
  }
  summaryVec <<- c(matched.len, perct, length(anot.id), sum(!hit.inx), ncol(dataSet$data.anot), ncol(dataSet$meta.info), sprintf("%4.2e", signif(totalCount ,3)), sprintf("%4.2e",signif(avgCount, 3)), sprintf("%4.2e",signif(minCount, 3)), sprintf("%4.2e",signif(maxCount,3)), lvls)  
  return(matched.len);   
}


PerformExpressNormalization <- function(norm.opt, var.thresh, count.thresh, abundance, filterUnmapped){
  
  print("normalizing ....");
  msg <- "Only features with annotations are kept for further analysis.";
  
  if(filterUnmapped == "false"){
    # need to update those with annotations
    data1 <- readRDS("data.proc.rds");
    anot.id <- readRDS("annotation.rds");;
    hit.inx <- !is.na(anot.id);
    rownames(data1)[hit.inx] <- anot.id[hit.inx];
    data1 <- RemoveDuplicates(data1, "mean", quiet=T);
    raw.data.anot <- data <- dataSet$data.anot <- data1;
  }else{
    raw.data.anot <- data <- readRDS("orig.data.anot");
  }
  
  if (dataSet$type == "count"){
    sum.counts <- apply(data, 1, sum, na.rm=TRUE);
    rm.inx <- sum.counts < count.thresh;
    data <- data[!rm.inx,];
    msg <- paste(msg, "Filtered ", sum(rm.inx), " genes with low counts.", collapse=" ");
  }else{
    avg.signal <- apply(data, 1, mean, na.rm=TRUE)
    p05 <- quantile(avg.signal, 0.05)
    all.rows = nrow(data)
    rm.inx = avg.signal < p05
    data <- data[!rm.inx,]
    msg <- paste(msg, "Filtered ", sum(rm.inx), " genes with low relative abundance (average expression signal).", collapse=" ");
  }
  
  data <- PerformDataNormalization(data, norm.opt);
  if(length(data)==1 && data == 0){
    return(0);
  }
  msg <- paste(norm.msg, msg);
  
  filter.val <- apply(data, 1, IQR, na.rm=T);
  nm <- "Interquantile Range";
  rk <- rank(-filter.val, ties.method='random');
  kp.pct <- (100 - var.thresh)/100;
  
  remain <- rk < nrow(data)*kp.pct;
  data <- data[remain,];
  msg <- paste(msg, paste("Filtered ", sum(!remain), " low variance genes based on IQR"), collapse=" ");
  
  dataSet$data.anot <- raw.data.anot[remain,]
  dataSet$data.norm <- data;
  
  # save normalized data for download user option
  write.csv(dataSet$data.norm, file="data_normalized.csv");
  
  current.msg <<- msg;
  dataSet <<- dataSet;
  
  processedObj <- list();#for omicsnet
  processedObj$name <- "rna_b_omicsanalyst.json"
  processedObj$type <- "rna.b"
  processedObj$data.proc <- dataSet$data.norm
  processedObj$feature.nms <- rownames(processedObj$data.proc)
  processedObj$sample.nms <- colnames(processedObj$data.proc)
  meta = dataSet$meta.info
  rownames(meta) <-  colnames(processedObj$data.proc)
  processedObj$meta <- meta
  library(RJSONIO)
  sink(processedObj$name);
  cat(toJSON(processedObj));
  sink();
  
  return(1);
}


# note, we do both filtering and normalization
PerformDataNormalization <- function(data, norm.opt){
  set.seed(1337);
  msg <- NULL;
  row.nms <- rownames(data);
  col.nms <- colnames(data);
  if(norm.opt=="log"){
    min.val <- min(data[data>0], na.rm=T)/10;
    numberOfNeg = sum(data<=0, na.rm = TRUE) + 1; 
    totalNumber = length(data)
    if((numberOfNeg/totalNumber)>0.2){
      msg <- paste(msg, "Can't perform log2 normalization, over 20% of data are negative. Try a different method or maybe the data already normalized?", collapse=" ");
      print(msg);
      norm.msg <<- current.msg <<- msg;
      return(0);
    }
    data[data<=0] <- min.val;
    data <- log2(data);
    msg <- paste(msg, "Log2 transformation.", collapse=" ");
  }else if(norm.opt=="vsn"){
    require(limma);
    data <- normalizeVSN(data);
    msg <- paste(msg, "VSN normalization.", collapse=" ");
  }else if(norm.opt=="quantile"){
    require('preprocessCore');
    data <- normalize.quantiles(data, copy=TRUE);
    msg <- paste(msg, "Quantile normalization.", collapse=" ");
  }else if(norm.opt=="combined"){
    require(limma);
    data <- normalizeVSN(data);
    require('preprocessCore');
    data <- normalize.quantiles(data, copy=TRUE);
    msg <- paste(msg, "VSN followed by quantile normalization.", collapse=" ");
  }else if(norm.opt=="logcount"){ # for count data, do it in DE analysis, as it is dependent on design matrix
    require(edgeR);
    nf <- calcNormFactors(data);
    y <- voom(data,plot=F,lib.size=colSums(data)*nf);
    data <- y$E; # copy per million
    msg <- paste(msg, "Limma based on log2-counts per million transformation.", collapse=" ");
  }else{
    # should do best guess for count data for plotting and filtering
    if(dataSet$type == "count"){
      if(sum(data > 100) > 100){ # now we think it is raw counts
        require(edgeR);
        nf <- calcNormFactors(data);
        y <- voom(data,plot=F,lib.size=colSums(data)*nf);
        data <- y$E; # copy per million
      }
    }
    msg <- paste(msg, "No log normalization was performed.", collapse=" ");
    print(msg);
  }
  norm.msg <<- msg;
  rownames(data) <- row.nms;
  colnames(data) <- col.nms;
  return(data);
}

# note, setup the main class, keep the original order
SetMainClass<-function(cls.lbl){
  lbls <- as.character(dataSet$meta.info[[cls.lbl]]);
  lvls.orig <- unique(lbls);
  cls <- factor(lbls, levels=lvls.orig, ordered=T);
  dataSet$cls <- cls; # record main cls
  dataSet <<- dataSet;
  return(levels(cls));
}

SetSelectedMetaInfo <- function(meta0, meta1, block1){
  if(meta0 == "NA"){
    return(0);
  }else{
    cls <- dataSet$meta.info[, meta0];
    dataSet$fst.cls <- cls; # for PCA plotting
    block <- NULL;
    dataSet$sec.cls <- "NA"
    if(meta1 != "NA"){
      if(block1){
        block <- dataSet$meta.info[, meta1]
      }else{ # two factor
        cls <- interaction(dataSet$meta.info[, c(meta0, meta1)], sep = ".", lex.order = TRUE);
      }
      dataSet$sec.cls <- dataSet$meta.info[, meta1]; # for pca coloring
    }
    dataSet$cls <- cls; # record main cls;
    dataSet$block <- block;
    dataSet <<- dataSet;
    return(levels(cls));
  }
}

SetupDesignMatrix<-function(deMethod){
  cls <- dataSet$cls; 
  design <- model.matrix(~ 0 + cls) # no intercept
  colnames(design) <- levels(cls);
  dataSet$design <- design;
  dataSet$de.method <- deMethod;
  dataSet <<- dataSet;
  
  return(1);
}

# perform differential analysis
# default: all pair-wise comparison (A-B) + (B-C) + (A-C)
# custom: only compare two groups (A-C)
# time: compare consecutive groups (B-A) + (C-B)
# reference: all others against common reference (A-C) + (B-C)
# nested: (A-B)+(C-D) 
PerformDEAnal<-function (anal.type = "default", par1 = NULL, par2 = NULL, nested.opt = "intonly"){
  set.seed(1337);
  par1 <<- par1
  par2 <<- par2
  nested.opt <<-nested.opt
  myargs <- list()
  cls <- dataSet$cls
  dataSet$comp.type = anal.type
  grp.nms <- levels(cls)
  
  if (anal.type == "default") {
    inx = 0
    for (m in 1:(length(grp.nms) - 1)) {
      for (n in (m + 1):length(grp.nms)) {
        inx <- inx + 1
        myargs[[inx]] <- paste(grp.nms[m], "-", grp.nms[n], sep = "");
      }
    }
    filename = "SigGene_pairwise";
  } else if (anal.type == "time") {
    for (i in 2:length(grp.nms)) {
      myargs[[i - 1]] <- paste(grp.nms[i], "-", grp.nms[i-1], sep = "")
    }
    filename = "SigGene_time_series"
  } else if (anal.type == "custom") {
    grp.nms <- strsplit(par1, " vs. ")[[1]]
    myargs[[1]] <- paste(grp.nms, collapse = "-")
    dataSet$grp.nms = grp.nms;
    filename = paste("SigGene_", paste(grp.nms, collapse = "_vs_"), sep = "")
  } else if (anal.type == "reference") {
    ref <- par1;
    cntr.cls <- grp.nms[grp.nms != ref]
    myargs <- as.list(paste(cntr.cls, "-", ref, sep = ""));
    filename = paste("SigGene_reference_", ref, sep = "");
  } else if (anal.type == "nested") {
    grp.nms1 <- strsplit(par1, " vs. ")[[1]]
    grp.nms2 <- strsplit(par2, " vs. ")[[1]]
    if (all(grp.nms1 == grp.nms2)) {
      current.msg <<- paste("The two nested groups are the same. Please choose two different groups.")
      return(0)
    }
    grp.nms <- unique(c(grp.nms1, grp.nms2))
    if (nested.opt == "intonly") {
      myargs[[1]] <- paste("(", paste(grp.nms1, collapse = "-"), ")-(", paste(grp.nms2, collapse = "-"), ")", sep = "")
    } else {
      myargs[[1]] <- paste(grp.nms1, collapse = "-")
      myargs[[2]] <- paste(grp.nms2, collapse = "-")
      myargs[[3]] <- paste("(", paste(grp.nms1, collapse = "-"), ")-(", paste(grp.nms2, collapse = "-"), ")", sep = "")
    }
    filename = paste("SigGene_nested_", paste(paste(grp.nms1, collapse = "_vs_"), "_", paste(grp.nms2, collapse = "_vs_"), sep = ""), sep = "")
  } else {
    print(paste("Not supported: ", anal.type))
  }
  
  library(limma)
  design <- dataSet$design
  myargs[["levels"]] <- design
  contrast.matrix <- do.call(makeContrasts, myargs)
  if (dataSet$de.method == "limma") {
    if (is.null(dataSet$block)) {
      fit = lmFit(dataSet$data.norm, design)
    } else {
      corfit <- duplicateCorrelation(dataSet$data.norm, design, block = dataSet$block)
      fit <- lmFit(dataSet$data.norm, design, block = dataSet$block, correlation = corfit$consensus)
    }
    
    if (!is.fullrank(design)) {
      current.msg <<- paste("This metadata combination is not full rank! Please use other combination.")
      return(0)
    }
    
    df.residual <- fit$df.residual
    if (all(df.residual == 0)) {
      current.msg <<- paste("There is not enough replicates in each group (no residual degrees of freedom)!")
      return(0);
    }
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    topFeatures <- topTable(fit2, number = Inf, adjust.method = "fdr");
    
  } else if (dataSet$de.method == "deseq2"){
    # only for small data set (< 60)
    if(length(cls) > 60){
      current.msg <<- "For large sample size (>60), use limma or edgeR."; 
      return(0);
      
    }else{ # use microservice
      # use microservice
      print("Peforming DESeq2 ....");
      
      library(RSclient);
      rsc <- RS.connect();
      RS.assign(rsc, "my.dir", getwd()); 
      RS.eval(rsc, setwd(my.dir));
      
      dat.out <- list(data=dataSet, contrast.matrix = contrast.matrix);
      RS.assign(rsc, "dat.in", dat.out); 
      my.fun <- function(){
        library(DESeq2);
        dataSet <- dat.in$data;
        contrast.matrix <- dat.in$contrast.matrix;
        
        if (is.null(dataSet$sec.cls)){
          colData <- data.frame(dataSet$fst.cls)
          colnames(colData) <- "condition"
          dds <- DESeqDataSetFromMatrix(countData=round(dataSet$data.anot), colData = colData, design = ~condition);
        } else {
          colData <- data.frame(dataSet$fst.cls, dataSet$sec.cls, dataSet$cls);
          colnames(colData) <- c("condition", "type", "condition_type");
          dds <- DESeqDataSetFromMatrix(countData=round(dataSet$data.anot), colData = colData, design = ~condition_type);
        }   
        
        dds <- DESeq(dds, betaPrior=TRUE)
        vec <- as.numeric(c(0, contrast.matrix[,1]));
        res <- results(dds, contrast = vec, independentFiltering = FALSE, cooksCutoff = Inf);
        topFeatures <- data.frame(res@listData);
        rownames(topFeatures) <- rownames(res);
        nms <- colnames(topFeatures);
        nms[which(nms == "padj")] <- "adj.P.Val";
        nms[which(nms == "pvalue")] <- "P.Value";
        nms[which(nms == "log2FoldChange")] <- "logFC";
        colnames(topFeatures) <- nms;
        topFeatures <- topFeatures[c(2,1,3,4,5,6)];
        # order the result based on raw p
        ord.inx <- order(topFeatures$P.Value);
        topFeatures <- topFeatures[ord.inx, ];
        return(topFeatures);
      }
      RS.assign(rsc, my.fun);
      topFeatures <-  RS.eval(rsc, my.fun());
      RS.close(rsc);
    }
  } else {
    library(edgeR)
    y <- DGEList(counts = dataSet$data.anot, group = dataSet$cls)
    y <- calcNormFactors(y)
    y <- estimateGLMCommonDisp(y, design, verbose = FALSE)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, contrast = contrast.matrix)
    topFeatures <- topTags(lrt, n = Inf)$table
    nms <- colnames(topFeatures)
    nms[which(nms == "FDR")] <- "adj.P.Val"
    colnames(topFeatures) <- nms
  }
  dataSet$filename <- filename;
  dataSet$resTable <- topFeatures;
  dataSet <<- dataSet;
  return(1)
}

# update result based on new cutoff
GetSigGenes <-function(res.nm, p.lvl, fc.lvl, update=T, inx){
  total = nrow(dataSet$resTable);
  resTable <- dataSet$resTable;
  filename <- dataSet$filename;
  filename <- paste(filename, "_", res.nm, ".csv", sep="");
  if (update){
    current.msg <<- "";
  }
  # select based on p-value
  if(dataSet$type == "array"){
    hit.inx.p <- resTable$adj.P.Val <= p.lvl; 
  } else {
    hit.inx.p <- resTable$adj.P.Val <= p.lvl; 
  }
  
  resTable<-resTable[hit.inx.p,,drop=F];
  if (nrow(resTable) == 0){
    current.msg <<- paste(current.msg, "No significant genes were identified using the given design and cutoff."); 
  }
  # now rank by logFC, note, the logFC for each comparisons 
  # are returned in resTable before the AveExpr columns 
  # for two-class, only one column, multiple columns can be involved
  # for > comparisons - in this case, use the largest logFC among all comparisons
  #if (fc.lvl > 0){ # further filter by logFC
  if (dataSet$de.method=="limma"){
    hit.inx <- which(colnames(resTable) == "AveExpr");
  } else if (dataSet$de.method=="deseq2"){
    hit.inx <- which(colnames(resTable) == "baseMean");  
  } else {
    hit.inx <- which(colnames(resTable) == "logCPM");
  }
  maxFC.inx <- hit.inx - 1; # not sure if this is also true for edgeR
  logfc.mat <- resTable[,1:maxFC.inx, drop=F];
  pos.mat <- abs(logfc.mat);
  fc.vec <- apply(pos.mat, 1, max);
  hit.inx.fc <- fc.vec >= fc.lvl;
  resTable<-resTable[hit.inx.fc,,drop=F];
  if (nrow(resTable) == 0){
    current.msg <<- paste(current.msg, "No significant genes were identified using the given design and cutoff."); 
    
  }
  #}
  
  ### Note, rowname of resTable must be entrez ID
  
  de.Num <- nrow(resTable);
  
  # display at most 5000 genes for the server (two main reasons)
  # 1) should not have more 22% (human: 23000) DE of all genes (biological)
  # 2) IE canvas can display no more than 6800 pixels (computational)
  if (nrow(resTable) > 5000){
    resTable <- resTable[1:5000,];
    current.msg <<- paste(current.msg, " Due to computational constraints, only top 5000 genes will be used. ", collapse="\n");
  }
  
  # may need to update data, class and meta.info
  data <- dataSet$data.norm;
  cls <- dataSet$cls; 
  meta.info <- dataSet$meta.info;
  grp.nms <- levels(cls);
  
  hit.inx <- cls %in% grp.nms;
  if (sum(hit.inx) < length(hit.inx)){
    current.msg <<- paste(current.msg, "Only groups selected for comparisons: ", paste(grp.nms, collapse=", "), "are included.");
    cls <- factor(cls[hit.inx]);
    cls.lvls <- levels(cls);
    data <- data[,hit.inx];
    meta.info <- dataSet$meta.info[hit.inx,];
  }
  
  saveRDS(data, file="data.stat");
  dataSet$resTable = dataSet$resTable[order(dataSet$resTable$adj.P.Val),] 
  dataSet$resTable = dataSet$resTable[which(!rownames(dataSet$resTable) %in% rownames(resTable)),]
  dataSet$resTable = rbind(resTable, dataSet$resTable);
  
  dataSet$sig.mat <- resTable;
  if (dataSet$annotated){ # annotated to entrez
    anot.id <- rownames(dataSet$resTable);
    gene.anot <- doEntrezIDAnot(anot.id);
    write.csv(cbind(EntrezID=anot.id, signif (dataSet$resTable,5), Symbols = gene.anot$symbol,  Name=gene.anot$name), row.names=F, file=filename);
  } else if (file.exists("annotation.rds")){ # annotation information available
    anot.id <- readRDS("annotation.rds");
    feature.vec <- rownames(dataSet$resTable);
    entrez.vec <- anot.id[feature.vec];
    gene.anot <- doEntrezIDAnot(entrez.vec);
    write.csv(cbind(signif (dataSet$resTable,5), EntrezID=entrez.vec, Symbols = gene.anot$symbol,  Name=gene.anot$name), row.names=F, file=filename);
    rownames(gene.anot) <- feature.vec;
  } else {
    gene.anot <- NULL;
    write.csv(signif(resTable,5), file=filename);
  }
  if(is.null(gene.anot)){
    dataSet$sig.genes.symbols <- rep("NA",nrow(resTable));
  }else{
    dataSet$sig.genes.symbols <- gene.anot$symbol;
  }
  dataSet$cls.stat <- cls;
  dataSet$meta.stat <- meta.info;
  
  # now do protein mapping for network only applicable for annotated
  
  dataSet$name <- res.nm;
  
  gene <- rownames(resTable);
  
  logFC <- unname(logfc.mat[,1]);
  geneList <- paste(gene, logFC, collapse="\n");
  up = nrow(resTable[which(logfc.mat[,selectedFactorInx]> fc.lvl),])
  down = nrow(resTable[which(logfc.mat[,selectedFactorInx]< -fc.lvl),])
  
  dataSet <<- dataSet;
  data.norm <- dataSet$data.norm
  colnames(data.norm) = NULL
  lst = list(colnames(dataSet$data.norm),data.norm, dataSet$meta.info, dataSet$resTable, rownames(data.norm), org=data.org)
  require(RJSONIO)
  json.obj <- toJSON(lst);
  sink("NetworkAnalyst_matrix.json");
  cat(json.obj);
  return(c(filename, de.Num, geneList, total, up, down));
}

GetExpressResultColNames<-function(){
  resT <- readRDS("ExpressResT.rda");
  colnames(resT);
}

GetExpressResultGeneIDs<-function(){
  return(rownames(dataSet$resTable));
}

GetExpressGeneIDType<-function(){
  return(dataSet$id.current);
}

GetExpressResultMatrix <-function(inxt){
  
  inxt = as.numeric(inxt)
  if (dataSet$de.method=="limma"){
    inx =match("AveExpr", colnames(dataSet$resTable))
  } else if (dataSet$de.method=="deseq2"){
    inx =match("baseMean", colnames(dataSet$resTable))
  } else {
    inx =match("logCPM", colnames(dataSet$resTable))
  }
  res = dataSet$resTable;
  res = res[,-(1:inx-1)]
  res <- cbind(dataSet$resTable[,inxt], res);
  colnames(res)[1] <- colnames(dataSet$resTable)[inxt];
  saveRDS(res, "ExpressResT.rda");
  return(signif(as.matrix(res), 5));
}


GetExpressResultGeneIDLinks <- function(){
  ids <- rownames(dataSet$resTable);
  symbs <- doEntrez2SymbolMapping(ids);
  annots <- paste("<a href='http://www.ncbi.nlm.nih.gov/gene?term=", ids, "' target='_blank'>", symbs, "</a>", sep="");
  return(annots);
}

GetExpressResultGeneSymbols<-function(){
  return(dataSet$sig.genes.symbols);
}

GetFactorNb<-function(){
  return(length(dataSet$meta.info));
}

PlotDataBox <- function(boxplotName, dpi, format){
  qc.boxplot(dataSet$data.norm, boxplotName, dpi, format);
}

PlotDataPCA <- function(pcaName, dpi, format,factor){
  qc.pcaplot(dataSet$data.norm, pcaName, dpi, format, factor);
}

PlotDataMeanStd <- function(densityName, dpi,format){
  qc.meanstd(dataSet$data.norm, densityName, dpi, format);
}


qc.density<- function(imgNm, dpi=72, format, factor){
  library("ggplot2")
  dat = dataSet$data.norm
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  dpi = as.numeric(dpi)
  
  df = data.frame(dataSet$data.norm, stringsAsFactors = FALSE)
  df = stack(df)
  sampleNms =gsub("-", ".", colnames(dataSet$data.norm))
  if(length(dataSet$meta.info) == 2){
    Factor1 = as.vector(dataSet$meta.info[,1])
    factorNm1 = colnames(dataSet$meta.info)[1]
    conv = data.frame(ind=sampleNms, Factor1=Factor1)
    colnames(conv) = c("ind", factorNm1);
    df1 = merge(df, conv, by="ind")
    Factor2 = as.vector(dataSet$meta.info[,2])
    factorNm2 = colnames(dataSet$meta.info)[2]
    conv = data.frame(ind=sampleNms, Factor2=Factor2)
    colnames(conv) = c("ind", factorNm2);
    df1 = merge(df1, conv, by="ind")
    df2 <- melt(df1, measure.vars=c(factorNm1,factorNm2))
    colnames(df2)[4] = "Conditions"
    g =ggplot(df2, aes(x=values)) + geom_line(aes(color=Conditions, group=ind), stat="density", alpha=0.6) + facet_grid(. ~ variable)
    width = 12
    height = 6
  }else{
    Conditions= as.character(dataSet$meta.info[,1]);
    conv = data.frame(ind=sampleNms, Conditions=Conditions)
    df1 = merge(df, conv, by="ind")
    g =ggplot(df1, aes(x=values)) + geom_line(aes(color=Conditions, group=ind), stat="density", alpha=0.6) 
    width = 8
    height = 6
  }
  Cairo(file=imgNm, width=width, height=height, type=format, bg="white", dpi=dpi, unit="in");
  print(g)
  dev.off();
}

qc.densitySample<- function(dat, imgNm, dpi=72, format){
  library("ggplot2")
  
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  dpi = as.numeric(dpi)
  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
  df = data.frame(dataSet$data.norm, stringsAsFactors = FALSE)
  df = stack(df)
  g = ggplot(df, aes(x=values, color=ind)) +
    geom_density()
  
  print(g)
  dev.off();
}



qc.meanstd <- function(dat, imgNm,dpi=72, format="png"){
  dpi = as.numeric(dpi)
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
  plot = meanSdPlot(dat, ranks=FALSE) 
  dev.off();
}

qc.boxplot <- function(dat, imgNm, dpi=72, format="png"){
  dpi = as.numeric(dpi)
  library('ggplot2')
  library('lattice');
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  subgene <- 10000;
  if (nrow(dat)>subgene) {
    set.seed(28051968);
    sg  <- sample(nrow(dat), subgene);
    Mss <- dat[sg,,drop=FALSE];
  } else {
    Mss <- dat;
  }
  
  subsmpl=100;
  if (ncol(Mss)>subsmpl) {
    set.seed(28051968);
    ss  <- sample(ncol(Mss), subsmpl)
    Mss <- Mss[,ss,drop=FALSE]
  } else {
    Mss <- Mss
  }
  
  sample_id <- rep(seq_len(ncol(Mss)), each = nrow(Mss));
  values  <- as.numeric(Mss)
  
  df = cbind(values, sample_id)
  
  df = data.frame(df)
  df$sample_id = factor(df$sample_id)
  xlower = unname(quantile(df$values, probs = c(0.01, 0.99), na.rm=TRUE)[1])
  xupper = unname(quantile(df$values, probs = c(0.01, 0.99), na.rm=TRUE)[2])
  height = length(unique(df$sample_id)) *20;
  if(height<450){
    height = 450
  }
  bp = ggplot(df, aes(sample_id, values)) +
    ylab("Values") + xlab("Samples") + scale_x_discrete(labels=colnames(dataSet$data.norm)) + ylim(xlower, xupper) + stat_boxplot(geom = "errorbar", color="black")+ geom_boxplot(outlier.size=0.5, outlier.alpha=0.4)
  bp = bp + coord_flip();
  
  Cairo(file=imgNm, width=600*dpi/72, height=height*dpi/72, unit="px",dpi=dpi, type=format, bg="white");
  print(bp);
  dev.off();
}

qc.pcaplot <- function(x, imgNm, dpi=72, format="png", factor){
  dpi = as.numeric(dpi);
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  require('lattice');
  require('ggplot2');
  require('reshape');
  pca <- prcomp(t(na.omit(x)));
  imp.pca<-summary(pca)$importance;
  xlabel = paste0("PC1"," (", 100*round(imp.pca[2,][1], 3), "%)")
  ylabel = paste0("PC2"," (", 100*round(imp.pca[2,][2], 3), "%)")
  names <- colnames(x);
  pca.res <- as.data.frame(pca$x);
  pca.res <- pca.res[,c(1,2)]
  # increase xlim ylim for text label
  xlim <- GetExtendRange(pca.res$PC1);
  ylim <- GetExtendRange(pca.res$PC2);
  if(length(dataSet$meta.info) == 2){
    Factor1 = as.vector(dataSet$meta.info[,1])
    factorNm1 = colnames(dataSet$meta.info)[1]
    pca.res[,factorNm1] = Factor1
    Factor2 = as.vector(dataSet$meta.info[,2])
    factorNm2 = colnames(dataSet$meta.info)[2]
    pca.res[,factorNm2] = Factor2
    pca.rest <- melt(pca.res, measure.vars=c(factorNm1,factorNm2))
    colnames(pca.rest)[4] = "Conditions"
    pca.rest$names = c(rownames(pca.res), rownames(pca.res))
    if(length(pca.rest$names)>20){
      pcafig = ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions, label=pca.rest$names)) +
        geom_point(size=3, alpha=0.5) + xlim(xlim)+ ylim(ylim) + xlab(xlabel) + ylab(ylabel) + facet_grid(. ~ variable)
    }else{
      require('ggrepel');
      pcafig = ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions, label=pca.rest$names)) +
        geom_point(size=4) + xlim(xlim)+ ylim(ylim) + xlab(xlabel) + ylab(ylabel) + geom_text_repel(force=1.5) + facet_grid(. ~ variable)
    }
    width = 12
    height = 6
  }else{
    Factor = dataSet$meta.info[,1];
    pca.rest = pca.res
    pca.rest$Conditions = Factor
    pca.rest$names = rownames(pca.res)
    if(length(rownames(pca.res))>20){
      pcafig = ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions)) +
        geom_point(size=3, alpha=0.5) + xlim(xlim)+ ylim(ylim)+ xlab(xlabel) + ylab(ylabel) 
    }else{
      require('ggrepel');
      pcafig = ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions, label=rownames(pca.res))) +
        geom_point(size=4) + xlim(xlim)+ ylim(ylim)+ xlab(xlabel) + ylab(ylabel) +geom_text_repel(force=1.5)+scale_color_manual(breaks=unique(pca.rest$Conditions), values=c("#00BFC4" ,"#F8766D"))
    }
    width = 8
    height = 6
  }
  Cairo(file=imgNm, width=width, height=height, type=format, bg="white", unit="in", dpi=dpi);
  print(pcafig);
  dev.off();
}


PlotLibSizeView<-function(imgNm,dpi=72, format="png",factor){
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  dpi=as.numeric(dpi)
  data_bef<-data.matrix(dataSet$data.anot);
  
  smpl.sums <- colSums(data_bef);
  
  library("ggplot2")
  data_bef<-data.matrix(dataSet$data.anot);
  smpl.sums <- colSums(data_bef);
  names(smpl.sums) <- colnames(data_bef);
  sampleNms = names(smpl.sums)
  df = data.frame(count=smpl.sums,ind=colnames(data_bef))
  
  if(length(dataSet$meta.info) == 2){
    Factor1 = as.vector(dataSet$meta.info[,1])
    factor1Nm = colnames(dataSet$meta.info)[1]
    conv = data.frame(ind=sampleNms, Factor1=Factor1)
    colnames(conv) = c("ind", factor1Nm)
    df1 = merge(df, conv, by="ind")
    Factor2 = as.vector(dataSet$meta.info[,2])
    factor2Nm = colnames(dataSet$meta.info)[2]
    conv = data.frame(ind=sampleNms, Factor2=Factor2)
    colnames(conv) = c("ind", factor2Nm)
    df1 = merge(df1, conv, by="ind")
    df2 <- melt(df1, measure.vars=c(factor1Nm,factor2Nm))
    colnames(df2)[4] = "Conditions"
    if(length(df2$ind)>20){
      g = ggplot(df2, aes(x = Conditions, y = count, fill=Conditions))+
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), dotsize=0.7) + ylab("Sum") + facet_grid(. ~ variable);
    }else{
      g = ggplot(df2, aes(x = Conditions, y = count, fill=Conditions, label=ind))+
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), dotsize=0.7) + geom_text_repel(force=5) + ylab("Sum") + facet_grid(. ~ variable);
    }
    width = 12
    height = 6
  }else{
    Conditions= as.character(dataSet$meta.info[,1]);
    conv = data.frame(ind=sampleNms, Conditions=Conditions)
    df1 = merge(df, conv, by="ind")
    if(length(df1$ind)>20){
      g = ggplot(df1, aes(x = Conditions, y = count, fill=Conditions))+
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), dotsize=0.7) + xlab("Sum");
    }else{
      g = ggplot(df1, aes(x = Conditions, y = count, label=ind, fill=Conditions))+
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), dotsize=0.7) + geom_text_repel(force=5) + xlab("Sum");
    }
    width = 8
    height = 6
  }
  
  Cairo(file=imgNm, width=width, height=height, unit="in", type=format, bg="white",dpi=dpi);
  
  #names(smpl.sums) <- colnames(data_bef);
  print(g);
  dev.off(); 
}

PlotMDS <- function(imgName, format){
  require(edgeR)
  require(RColorBrewer)
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  Cairo(file=imgName, width=580, type=format, bg="white",dpi=72)
  levels(dataSet$cls) <- brewer.pal(nlevels(dataSet$cls), "Set1")
  col.group <- dataSet$cls
  col.group <- as.character(col.group)
  plotMDS(dataSet$data.norm, col=col.group, xlab = "Dimension 2", ylab = "Dimension 1")
  title(main="MDS")
  dev.off();
} 

GetSummaryData <- function(){
  return(summaryVec);
}

GetDensityPlot <- function() {
  data = list();
  dat = dataSet$data.norm;
  for(i in 1: ncol(dat)){
    data[[i]] = dat[,i]
    names(data)[i] = colnames(dat)[i];
  }
  densityList = list();
  for(i in 1: length(data)){
    d = density(data[[i]])
    df = data.frame(d$x, d$y)
    colnames(df) = c("x","y");
    densityList[[i]] = df
  }
  names(densityList) = colnames(dataSet$data.norm);
  jsonNm <- "density.json";
  lst = list()
  
  lst=list(
    density= densityList,
    class= dataSet$cls
  )
  require(RJSONIO);
  json.obj <- toJSON(lst);
  sink(jsonNm);
  cat(json.obj);
}

PlotMAPlot<- function(imgNm, dpi=72, format,pvalue, fc, inx){
  library('ggplot2');
  dpi=as.numeric(dpi);
  inx = as.numeric(inx);
  pvalue=as.numeric(pvalue);
  fc=as.numeric(fc);
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  Cairo(file=imgNm, width=700, height=560,type=format, bg="white",dpi=100)
  pvalue = as.numeric(pvalue);
  fc = as.numeric(fc);
  res = dataSet$resTable
  if(dataSet$de.method =="deseq2"){
    #res = res[!apply(sapply(res, function(x) abs(scale(x)) >= 5), 1, any), ];
    res['log2(baseMean)'] = log2(res$baseMean);
  }
  
  # select based on p-value
  if(dataSet$type == "array"){
    res$significant = ifelse(res[,inx] > fc & res$P.Value < pvalue , 1, ifelse(res[,inx] < -fc & res$P.Value < pvalue, -1, 0))
  } else {
    res$significant = ifelse(res[,inx] > fc & res$adj.P.Val < pvalue , 1, ifelse(res[,inx] < -fc & res$adj.P.Val < pvalue, -1, 0))
  }
  res$significant = as.factor(res$significant)
  yCol= colnames(res)[inx]
  
  if (dataSet$de.method=="limma"){
    maplot = ggplot(res, aes_string(x="AveExpr", y=yCol, color="significant"))
  } else if (dataSet$de.method=="deseq2"){
    maplot = ggplot(res, aes_string(x="log2(baseMean)", y=yCol, color="significant"))
  } else {
    maplot = ggplot(res, aes_string(x="logCPM", y=yCol, color="significant"))
  }
  
  maplot = maplot +
    geom_point(size=0.5, alpha=0.5) +
    geom_hline(color = "blue3", yintercept = 0) +
    stat_smooth(se = FALSE, method = "loess", color = "red3") +
    scale_color_manual(values=c("-1" = "green", "0" = "black", "1" = "red"))+ theme(legend.position="none")
  print(maplot)
  dev.off()
}

GetMetaColLength<- function(){
  if (dataSet$de.method=="limma"){
    inx =match("AveExpr", colnames(dataSet$resTable))
  } else if (dataSet$de.method=="deseq2"){
    inx =match("baseMean", colnames(dataSet$resTable))
  } else {
    inx =match("logCPM", colnames(dataSet$resTable))
  }
  resT = dataSet$resTable;
  resT = resT[,1:inx-1]
  return(length(colnames(resT)));
}


GetMetaInfoLength<- function(){
  return(length(dataSet$meta.info));
}

GetMetaCol<- function(){
  colNms <- colnames(dataSet$resTable);
  if (dataSet$de.method=="limma"){
    inx =match("AveExpr", colNms)
  } else if (dataSet$de.method=="deseq2"){
    inx =match("baseMean", colNms)
  } else {
    inx =match("logCPM", colNms)
  }
  resT <- dataSet$resTable;
  resT <- resT[,1:inx-1]
  return(colnames(resT));
}

Volcano.Anal <- function(paired=FALSE, fcthresh, threshp, analType, inx){
  inx = as.numeric(inx)
  print("Prepare volcano anal");
  if(anal.type == "metadata"){
    if(dataSet$name != selDataNm){
      dataSet <- readRDS(selDataNm);
    }
    data <- as.matrix(inmex.ind[selDataNm][[1]])
    p.value <- data[, "Pval"]
    fcthresh = 0;
    
  }else{
    data <- as.matrix(dataSet$resTable);
    
    if(dataSet$type == "array"){
      p.value <- data[, "adj.P.Val"];
    } else {
      p.value <- data[, "adj.P.Val"];
    }
  }
  fcthreshu <<- fcthresh
  
  if (analType == "qPCR"){
    inx.p <- p.value < 1;
  } else {
    inx.p <- p.value <= threshp;
  }
  zero.inx <- p.value == 0;
  if(sum(zero.inx)>0){
    p.value[zero.inx] <- min(p.value[!zero.inx])/10;
  }
  p.log <- -log10(p.value);
  
  if (dataSet$annotated){ # annotated to entrez
    anot.id <- rownames(data);
    gene.anot <- doEntrezIDAnot(anot.id);
  }else{
    anot.id <- rownames(data);
    gene.anot <- data.frame(gene_id=anot.id, symbol=anot.id, stringsAsFactors=FALSE)
    init.lib <<- "NA"
  }
  
  #gene symbol to be used for boxplot   
  
  # create a named matrix of sig vars for display
  fc.log <- data[, inx];
  hit.maxPos <- (which(fc.log> 10) )
  hit.maxNeg <- (which(fc.log< -10) )
  fc.log[hit.maxPos] = 10;
  fc.log[hit.maxNeg] = 10;
  #fc.all <- res$fc.all;
  
  if(fcthresh != 0){
    inx.up = fc.log > fcthresh & p.value < threshp;
    inx.down = fc.log < -fcthresh & p.value < threshp;
  }else{
    inx.up = fc.log > 0 & p.value < threshp;
    inx.down = fc.log < 0 & p.value < threshp;
  }
  
  # create named sig table for display
  inx.imp <- (inx.up | inx.down) & inx.p;
  sig.var <- cbind(fc.log[inx.imp,drop=F], p.value[inx.imp,drop=F], p.log[inx.imp,drop=F]);
  colnames(sig.var) <- c("log2(FC)", "p.value", "-log10(p)");
  # first order by log(p), then by log(FC)
  ord.inx <- order(sig.var[,3], abs(sig.var[,1]), decreasing=T);
  sig.var <- sig.var[ord.inx,, drop=F];
  
  sig.var <- signif (sig.var, 5);
  sig.var1 <- sig.var;
  sig.var1 = cbind(rownames(sig.var), sig.var);
  colnames(sig.var1) <- c("name", "log2(FC)", "p.value", "-log10(p)");
  
  ###########################
  ## for Volcano data
  ##########################
  
  if(init.lib != "NA"){
    PerformVolcanoEnrichment("abc", init.lib, "null", "all", inx)
  }
  
  fileName <- "volcano.csv";
  jsonNm <- "volcano.json";
  require(RJSONIO);
  json.obj <- toJSON(sig.var1);
  sink(jsonNm);
  cat(json.obj);
  sink();
  write.csv(signif (sig.var,5),file=fileName);
  colnames(gene.anot)[1] = "anot.id"
  volcano <- list (
    raw.threshx = fcthresh,
    raw.threshy = threshp,
    paired = paired,
    thresh.y = -log10(threshp),
    fc.symb =rownames(data),
    fc.log = fc.log,
    fc.log.uniq = jitter(fc.log),
    inx.up = inx.up,
    inx.down = inx.down,
    p.log = p.log,
    inx.p = inx.p,
    sig.mat = sig.var,
    conv = gene.anot
  );
  
  require(RJSONIO);
  json.obj <- toJSON(volcano);
  sink("volcano2.json");
  cat(json.obj);
  sink();
  
  require(RJSONIO);
  
  if(init.lib == "NA"){
    enr.mat = "NA"
  }
  write.csv(enr.mat, file="enrichment_result.csv", row.names=T);
  sink("enrichment_result.json");
  cat(json.obj);
  sink();
  
}

# perform limma on given two groups selected 
# used by integarative analysis
PerformLimmaDE<-function(dataName, grps, p.lvl, fc.lvl=NULL){
  
  print("doing differential analysis ....");
  dataSet <- readRDS(dataName);
  dataSet$pval <- p.lvl
  if(length(levels(dataSet$cls))>2){ 
    grp.nms <- strsplit(grps, " vs. ")[[1]];
    sel.inx <- as.character(dataSet$cls) %in% grp.nms;
  }else{
    sel.inx <- rep(T, ncol(dataSet$data));
  }
  
  group <- factor(dataSet$cls[sel.inx]); # note regenerate factor to drop levels 
  data <- dataSet$data[, sel.inx];
  
  res.limma <- PerformLimma(data, group);
  res.all <- GetLimmaResTable(res.limma$fit.obj);
  
  if(!is.null(fc.lvl)){
    hit.inx <- abs(res.all$logFC)>= fc.lvl & res.all$adj.P.Val <= p.lvl
  }else{
    hit.inx <- res.all$adj.P.Val <= p.lvl
  }
  # note, hit.inx can contain NA, not T/F
  hit.inx <- which(hit.inx);
  res<-res.all[hit.inx,];
  
  # rm .txt suffix for new names
  shortNm <- substring(dataName, 0, nchar(dataName)-4);
  write.csv(signif(res[,-1],5), file=paste("SigGenes_", shortNm, ".csv",sep=""));
  
  sig.count <- nrow(res);
  de.genes <- rownames(res);
  res.mat <- cbind(res.all$logFC, res.all$adj.P.Val);
  rownames(res.mat) <- rownames(res.all);
  non.sig.count <- nrow(data)-sig.count;
  rm(res.all);
  
  gc();
  RegisterData(dataSet);
  # record the sig gene vec
  return (c(1, sig.count, non.sig.count));
}

# perfor differential analysis for array/RNA seq data
# for two groups only (used for meta-analysis)
PerformLimma<-function(data, group){
  require(limma);
  data <- data;
  design <- model.matrix(~-1 + group);
  fit = lmFit(data, design)
  
  grps.cmp <- paste("group", levels(group)[2], " - ", "group", levels(group)[1], sep="");
  myargs <- list(grps.cmp, levels = design);
  contrast.matrix <- do.call(makeContrasts, myargs);
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit);
  gc();
  return (list(fit.obj=fit));
}

# get result table from eBayes fit object
GetLimmaResTable<-function(fit.obj){
  resTable <- topTable(fit.obj, number=Inf, adjust.method="BH");
  if(!is.null(resTable$ID)){ # for older version
    rownames(resTable) <- resTable$ID;
    resTable$ID <- NULL;
  }
  return (resTable);
}

# given a gene id, plot its expression profile as box plot
PlotSelectedGene<-function(gene.id, type){
  
  imgName <- paste("Gene_", gene.id, ".png", sep="");
  require(lattice);
  if(anal.type == "onedata"){
    ids <- rownames(dataSet$resTable);
    inx <- which(ids == gene.id);
    symb <- dataSet$sig.genes.symbols[inx]; 
    if(type== "volcano"){
      symb = "";
    }    
    if(dataSet$comp.type == "custom"){
      Cairo(file = imgName, width=280, height=320, type="png", bg="white");
      grp.nms = dataSet$grp.nms;
      inx = dataSet$cls %in% grp.nms;
      cls = dataSet$cls[inx]
      dat = dataSet$data.norm[,inx];
      myplot <- bwplot(dat[gene.id,] ~ as.character(cls), fill="#0000ff22", scales=list(x=list(rot=30)),
                       xlab="Class", ylab="Expression Pattern", main=symb);
    }else if(length(dataSet$sec.cls)>1){
      out.fac <- as.character(dataSet$sec.cls)
      in.fac <- as.character(dataSet$fst.cls)
      xlab = colnames(dataSet$meta.info[,1]);
      
      Cairo(file = imgName, dpi=72, width=320, height=320, type="png", bg="white");
      #ylim.ext <- GetExtendRange(dataSet$data.norm[gene.id, ], 12);
      layout <- c(2, 1);
      myplot<- bwplot(dataSet$data.norm[gene.id, ] ~ in.fac | out.fac, 
                      xlab="Factors", ylab="Expression Pattern", main=symb, scales=list(x=list(rot=30)),
                      fill="#0000ff22", layout=layout);
    }else{
      Cairo(file = imgName, width=280, height=320, type="png", bg="white");
      
      myplot <- bwplot(dataSet$data.norm[gene.id,] ~ as.character(dataSet$cls), fill="#0000ff22", scales=list(x=list(rot=30)),
                       xlab="Class", ylab="Expression Pattern", main=symb);
    }
    
  }else{ # metadata
    
    inmex.meta <- readRDS("inmex_meta.rds");
    if(inmex.meta$id.type == "entrez"){
      symb <- inmex.meta$gene.symbls[gene.id];
    }else{
      symb <- gene.id;
    }
    num <- sum(mdata.all == 1);
    # calculate width based on the dateset number
    if(num == 1){
      Cairo(file = imgName, width=280, height=320, type="png", bg="white");
      myplot <- bwplot(inmex.meta$plot.data[gene.id,] ~ as.character(inmex.meta$cls.lbl), fill="#0000ff22",
                       xlab="Class", ylab="Expression Pattern", main=symb, scales=list(x=list(rot=30)))
    }else{
      # calculate layout
      if(num < 6){
        layout <- c(num, 1);
        height=320;
        width=160*num;
      }else{
        rn <- round(num/2);
        layout <- c(rn, 2);
        height=500;
        width=160*rn;
      }
      
      Cairo(file = imgName, width=width, height=height, type="png", bg="white");
      data.lbl <- as.character(inmex.meta$data.lbl);
      data.lbl <- substr(data.lbl, 0, nchar(data.lbl)-4);
      
      # get counts in each data, same order as a levels
      counts <- table(data.lbl);
      # back to factor 
      data.lbl <- factor(data.lbl);
      
      # get new lbls to cut potential long names, and add sample numbers
      nlbls <- data.lbl;
      levels(nlbls) <- abbreviate(levels(nlbls),9);
      nlbls <- paste(levels(nlbls), "( n=", as.vector(counts), ")");
      # update labels
      data.lbl <- factor(data.lbl, labels=nlbls);
      # some time the transformed plot.data can switch class label, use the original data, need to be similar scale
      myplot <- bwplot(inmex.meta$plot.data[gene.id,] ~ as.character(inmex.meta$cls.lbl) | data.lbl, 
                       xlab="Datasets", ylab="Expression Pattern", main=symb, scales=list(x=list(rot=30)),
                       fill="#0000ff22", layout=layout);
    }
  }
  print(myplot); 
  dev.off();
}

PlotSelectedGeneLoading<-function(gene.id){
  if(anal.type == "metadata"){
    PlotSelectedGeneMeta(gene.id);
  }else{
    PlotSelectedGene(gene.id, "notVolcano");
  }
}

PlotSelectedGeneMeta<-function(gene.id){
  
  # first get gene symbol
  inmex.meta <- readRDS("inmex_meta.rds");
  if(inmex.meta$id.type == "entrez"){
    symb <- inmex.meta$gene.symbls[gene.id];
  }else{
    symb <- gene.id;
  }
  
  imgName <- paste("Gene_", gene.id, ".png", sep="");
  require(lattice);
  
  num <- sum(mdata.all == 1);
  # calculate width based on the dateset number
  if(num == 1){
    Cairo(file = imgName, width=280, height=320, type="png", bg="white");
    myplot <- bwplot(inmex.meta$plot.data[gene.id,] ~ as.character(inmex.meta$cls.lbl), fill="#0000ff22",
                     xlab="Class", ylab="Expression Pattern", main=symb, scales=list(x=list(rot=30)))
  }else{
    # this is a single long list 
    layout <- c(1, num);
    height <- 200*num;
    width <- 280;
    
    Cairo(file = imgName, width=width, height=height, type="png", bg="white");
    data.lbl <- as.character(inmex.meta$data.lbl);
    data.lbl <- substr(data.lbl, 0, nchar(data.lbl)-4);
    
    # get counts in each data, same order as a levels
    counts <- table(data.lbl);
    # back to factor 
    data.lbl <- factor(data.lbl);
    
    # get new lbls to cut potential long names, and add sample numbers
    nlbls <- data.lbl;
    levels(nlbls) <- abbreviate(levels(nlbls),9);
    nlbls <- paste(levels(nlbls), "( n=", as.vector(counts), ")");
    # update labels
    data.lbl <- factor(data.lbl, labels=nlbls);
    # some time the transformed plot.data can switch class label, use the original data, need to be similar scale
    myplot <- bwplot(inmex.meta$plot.data[gene.id,] ~ as.character(inmex.meta$cls.lbl) | data.lbl, 
                     xlab="Datasets", ylab="Expression Pattern", main=symb, scales=list(x=list(rot=30)),
                     fill="#0000ff22", layout=layout);
  }
  
  print(myplot); 
  dev.off();
}

PlotCmpdView <-function(cmpdNm, format="png", dpi=72, width=NA){
  if(anal.type == "onedata"){
    datanorm = dataSet$data.norm
  }else{
    datanorm = dataSet$data
  }
  clslbl = dataSet$meta.info[,1];
  imgName <- gsub("\\/", "_",  cmpdNm);
  imgName <- paste(imgName, "_dpi", dpi, ".", format, sep="");
  #indx<-which(rownames(boxplot_id)==cmpdNm);
  #gene.id <- boxplot_id[indx,1];
  gene.symb <<- doEntrez2SymbolMapping(cmpdNm);
  Cairo(file = imgName, dpi=dpi, width=230, height=230, type=format, bg="transparent");
  par(mar=c(4,3,1,2), oma=c(0,0,1,0));
  boxplot(datanorm[which(rownames(datanorm)==as.character(cmpdNm)),]~clslbl,las=2,col= unique(GetColorSchema(clslbl)));
  title(main=gene.symb, out=T);
  dev.off();
  return(imgName);
}

# retrun the json obj
SaveHeatmapJSON <- function(fileName){
  
  if(anal.type == "metadata"){
    json.res <- PrepareMetaHeatmapJSON();
  }else{
    json.res <- PrepareExpressHeatmapJSON();
  }
  
  require(RJSONIO);
  json.mat <- toJSON(json.res, .na='null');
  sink(fileName);
  cat(json.mat);
  sink();
  current.msg <<- "Data is now ready for heatmap visualization!";
  return(1);
}


# prepare data for heatmap plotting include 
# 1. meta info (top bars)
# 2. expression matrix
# 3. function annotation (left bars)
# all matrix will be combined, 
# 1 and 2 separated by a row of 'null' 
# 3 and 1+2 separated by a column of 'null'
PrepareExpressHeatmapJSON <- function(){
  sig.ids <- rownames(dataSet$sig.mat);
  stat.pvals <- dataSet$sig.mat$adj.P.Val;
  stat.fc <- dataSet$sig.mat$logFC; 
  
  # scale each gene 
  data.stat <- readRDS("data.stat");
  hit.inz = sig.ids %in% rownames(data.stat);
  sig.ids = sig.ids[hit.inz];
  dat <- t(scale(t(data.stat[sig.ids, , drop=F])));
  
  # now pearson and euclidean will be the same after scaleing
  dat.dist <- dist(dat); 
  
  orig.smpl.nms <- colnames(dat);
  orig.gene.nms <- rownames(dat);
  
  # do clustering and save cluster info
  # convert order to rank (score that can used to sort) 
  if(nrow(dat)> 1){
    dat.dist <- dist(dat);
    gene.ward.ord <- hclust(dat.dist, "ward.D")$order;
    gene.ward.rk <- match(orig.gene.nms, orig.gene.nms[gene.ward.ord]);
    gene.ave.ord <- hclust(dat.dist, "ave")$order;
    gene.ave.rk <- match(orig.gene.nms, orig.gene.nms[gene.ave.ord]);
    gene.single.ord <- hclust(dat.dist, "single")$order;
    gene.single.rk <- match(orig.gene.nms, orig.gene.nms[gene.single.ord]);
    gene.complete.ord <- hclust(dat.dist, "complete")$order;
    gene.complete.rk <- match(orig.gene.nms, orig.gene.nms[gene.complete.ord]);
    
    dat.dist <- dist(t(dat));
    smpl.ward.ord <- hclust(dat.dist, "ward.D")$order;
    smpl.ward.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ward.ord])
    smpl.ave.ord <- hclust(dat.dist, "ave")$order;
    smpl.ave.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ave.ord])
    smpl.single.ord <- hclust(dat.dist, "single")$order;
    smpl.single.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.single.ord])
    smpl.complete.ord <- hclust(dat.dist, "complete")$order;
    smpl.complete.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.complete.ord])
  }else{
    # force not to be single element vector which will be scaler
    stat.pvals <- matrix(stat.pvals);
    stat.fc <- matrix(stat.fc);
    gene.ward.rk <- gene.ave.rk <- gene.single.rk <- gene.complete.rk <- matrix(1);
    smpl.ward.rk <- smpl.ave.rk <- smpl.single.rk <- smpl.complete.rk <- 1:ncol(dat);
  }
  
  gene.cluster <- list(
    pval = stat.pvals, 
    fc = stat.fc,
    ward = gene.ward.rk,
    average = gene.ave.rk,
    single = gene.single.rk,
    complete = gene.complete.rk
  );
  
  sample.cluster <- list(
    ward = smpl.ward.rk,
    average = smpl.ave.rk,
    single = smpl.single.rk,
    complete = smpl.complete.rk
  );
  
  # prepare meta info    
  # 1) convert meta.data info numbers
  # 2) match number to string (factor level)
  meta <- data.frame(dataSet$meta.stat);
  grps <- colnames(meta);
  nmeta <- meta.vec <- NULL;
  uniq.num <- 0;
  for (i in 1:ncol(meta)){
    cls <- meta[,i];
    grp.nm <- grps[i];
    meta.vec <- c(meta.vec, as.character(cls))
    # make sure each label are unqiue across multiple meta data
    ncls <- paste(grp.nm, as.numeric(cls)); # note, here to retain ordered factor
    nmeta <- c(nmeta, ncls);
  }
  
  # convert back to numeric 
  nmeta <- as.numeric(as.factor(nmeta))+99;
  unik.inx <- !duplicated(nmeta)   
  
  # get corresponding names
  meta_anot <- meta.vec[unik.inx]; 
  names(meta_anot) <- nmeta[unik.inx]; # name annotatation by their numbers
  
  nmeta <- matrix(nmeta, ncol=ncol(meta), byrow=F);
  colnames(nmeta) <- grps;
  
  # for each gene/row, first normalize and then tranform real values to 30 breaks 
  res <- t(apply(dat, 1, function(x){as.numeric(cut(x, breaks=30))}));
  
  # note, use {} will lose order; use [[],[]] to retain the order
  
  
  if(dataSet$annotated){
    anot.id <- rownames(res);
    anot.res <- doEntrezIDAnot(anot.id);
    # single element vector will be converted to scalar, not array, need to prevent that
    gene.id = anot.res$symbol; if(length(gene.id) ==1) { gene.id <- matrix(gene.id) };
    gene.entrez = anot.res$gene_id; if(length(gene.entrez) ==1) { gene.entrez <- matrix(gene.entrez) };        
    gene.name = anot.res$name; if(length(gene.name) ==1) { gene.name <- matrix(gene.name) };
    
    json.res <- list(
      data.type = dataSet$type, 
      gene.id = anot.res$symbol,
      gene.entrez = gene.entrez,
      gene.name = anot.res$name,
      gene.cluster = gene.cluster,
      sample.cluster = sample.cluster,
      sample.names = orig.smpl.nms,
      meta = data.frame(nmeta),
      meta.anot = meta_anot,
      data = res
    );
  }else if(file.exists("annotation.rds")){
    # special gene.id and new gene.symbol
    anot.id <- rownames(res);
    anot.res <- doEntrezIDAnot(anot.id);
    gene.id = rownames(anot.res); if(length(gene.id) ==1) { gene.id <- matrix(gene.id) };
    gene.entrez = anot.res$gene_id; if(length(gene.entrez) ==1) { gene.entrez <- matrix(gene.entrez) };  
    gene.name = paste(anot.res$symbol, anot.res$name, sep=" | "); if(length(gene.name) ==1) { gene.name <- matrix(gene.name) };
    
    json.res <- list(
      data.type = dataSet$type, 
      gene.id = gene.id,
      gene.entrez = gene.entrez,
      gene.name = gene.name,
      gene.cluster = gene.cluster,
      sample.cluster = sample.cluster,
      sample.names = orig.smpl.nms,
      meta = data.frame(nmeta),
      meta.anot = meta_anot,
      data = res
    );
  }else{          
    gene.id = orig.gene.nms; if(length(gene.id) ==1) { gene.id <- matrix(gene.id) };
    json.res <- list(
      data.type = dataSet$type,
      gene.id = gene.id,
      gene.entrez = gene.id,
      gene.name = gene.id,
      gene.cluster = gene.cluster,
      sample.cluster = sample.cluster,
      sample.names = orig.smpl.nms,
      meta = data.frame(nmeta),
      meta.anot = meta_anot,
      data = res
    );
  }
  return(json.res);
}

SaveClusterJSON <- function(fileNm, clustOpt, opt){
  if(anal.type == "onedata"){
    SaveExpressClusterJSON(fileNm, clustOpt,opt);
  }else{
    initmetaloading <<- TRUE;
    SaveMetaClusterJSON(fileNm, clustOpt, opt);
  }
}

SaveClusterJSONLoading <- function(fileNm, clustOpt, nb){
  if(anal.type == "onedata"){
    SaveExpressClusterLoadingJSON(fileNm, clustOpt, nb);
  }else{
    SaveMetaClusterLoadingJSON(fileNm, clustOpt, nb);
  }
}

SaveExpressClusterLoadingJSON <- function(fileName, clustOpt, nb){
  dat <- dataSet$data.norm;
  pca3d <- list();
  dat <- na.omit(dat);
  nb = as.numeric(nb)
  if(clustOpt == "pca"){
    pca <- prcomp(t(dat), center=T, scale=T);
    imp.pca<-summary(pca)$importance;
    pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
    coords <- data.frame(t(signif(pca$rotation[,1:3], 5)));
    
    colnames(coords) <- NULL; 
    pca3d$score$xyz <- coords;
    pca3d$score$name <- doEntrez2SymbolMapping(rownames(pca$rotation));
    pca3d$score$entrez <-rownames(pca$rotation);
    weights = imp.pca[2,][1:3]
    mypos <- t(coords);
    meanpos = apply(abs(mypos),1, function(x){weighted.mean(x, weights)})
    df = data.frame(pos = meanpos, inx = seq.int(1,length(meanpos)))
    df = df[order(-df$pos),]
    if(nrow(df) >2000){
      inx = df$inx[c(1:nb)]
      mypos = mypos[inx,];
      pca3d$score$xyz = coords[inx]
      pca3d$score$name = pca3d$score$name[inx]
      pca3d$score$entrez = pca3d$score$entrez[inx]
    }
  }
  
  pca3d$cls = dataSet$meta.info;
  colnames(mypos) <- paste("Dim", 1:3, sep="");
  # see if there is secondary
  loadEntrez <<- pca3d$score$entrez
  rownames(mypos) = pca3d$score$name;
  
  write.csv(mypos, file="networkanalyst_3d_load_pos.csv");
  require(RJSONIO);
  json.mat <- toJSON(pca3d, .na='null');
  sink(fileName);
  cat(json.mat);
  sink();
  current.msg <<- "Annotated data is now ready for PCA 3D visualization!";
  return(1);
}

# single expression data
SaveExpressClusterJSON <- function(fileName, clustOpt, opt){
  dat <- dataSet$data.norm;
  pca3d <- list();
  dat <- na.omit(dat);
  
  if(clustOpt == "pca"){
    if(opt == "all"){
      pca <- prcomp(t(dat), center=T, scale=T);
    }else{
      dat = dat[which(rownames(dat) %in% loadEntrez),]
      pca <- prcomp(t(dat), center=T, scale=T);
    }
    imp.pca<-summary(pca)$importance;
    pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
    coords <- data.frame(t(signif(pca$x[,1:3], 5)));
  }else{ # tsne
    require('Rtsne');
    dat <- as.matrix(t(dat));
    max.perx <- floor((nrow(dat)-1)/3);
    if(max.perx > 30){
      max.perx <- 30;
    }
    res <- Rtsne(dat, dims = 3, perplexity=max.perx);
    pca3d$score$axis <- paste("t-SNE dim ", 1:3, sep="");
    coords <- data.frame(t(signif(res$Y, 5)));
  }
  
  colnames(coords) <- NULL; 
  pca3d$score$xyz <- coords;
  pca3d$score$name <- colnames(dataSet$data.norm);
  
  facA <- as.character(dataSet$fst.cls);
  if(all.numeric(facA)){
    facA <- paste("Group", facA);
  }
  pca3d$score$facA <- facA;
  
  mypos <- t(coords);
  colnames(mypos) <- paste("Dim", 1:3, sep="");
  # see if there is secondary
  if(length(dataSet$sec.cls) > 1){
    facB <- as.character(dataSet$sec.cls);
    if(all.numeric(facB)){
      facB <- paste("Group", facB);
    }
    pca3d$score$facB <- facB;
    
    # set shape based on the first group
    pca3d$score$shapes <- c("sphere", "triangle");
    
    # now set color based on 2nd group
    cols <- unique(GetColorSchema(dataSet$sec.cls));
    rgbcols <- col2rgb(cols);
    cols <- apply(rgbcols, 2, function(x){paste("rgb(", paste(x, collapse=","), ")", sep="")});
    pca3d$score$colors <- cols;
    
    mypos <- data.frame(factorA=facA, factorB=facB, mypos);
  }else{
    # now set color based on first group
    cols <- unique(GetColorSchema(dataSet$fst.cls));
    rgbcols <- col2rgb(cols);
    cols <- apply(rgbcols, 2, function(x){paste("rgba(", paste(x, collapse=","), ",1)", sep="")});
    pca3d$score$colors <- cols;
    mypos <- data.frame(factorA=facA, mypos);
  }
  pca3d$cls = dataSet$meta.info;
  rownames(mypos) = colnames(dataSet$data.norm);
  
  write.csv(mypos, file="networkanalyst_3d_pos.csv");
  require(RJSONIO);
  json.mat <- toJSON(pca3d, .na='null');
  sink(fileName);
  cat(json.mat);
  sink();
  current.msg <<- "Annotated data is now ready for PCA 3D visualization!";
  return(1);
}

SetVolcanoHigh<- function(ids){
  idsu <<- ids
  gene.vec <- unlist(strsplit(ids, "; "));
  volcanoHlVec <<- gene.vec;
  if(length(volcanoHlVec)>0){
    return(1);
  }else{
    return(0);
  }
}

# prepare seeds from metaanalysis result
.prepareExpressSeeds <- function(){
  
  gene.list <- list();
  if(anal.type == "metadata"){
    if(selectedNetDataset == "meta_dat"){
      gene.mat <- meta.mat;
      if(inmex.method != "votecount"){
        gene.list$metadata <- list(gene=rownames(gene.mat),logFC=unname(gene.mat[,1]), adjP = unname(gene.mat[,2]));   
      }else{
        gene.list$metadata <- list(gene=rownames(gene.mat),logFC=unname(gene.mat[,1]), adjP = unname(gene.mat[,1]));   
      }
    }else{
      #dataSet <- readRDS(selectedNetDataset);
      fit.obj.nm <- paste(selectedNetDataset, "fit.obj", sep=".");
      fit2i <- readRDS(fit.obj.nm);
      gene.mat <- GetLimmaResTable(fit2i);
      gene.mat <- gene.mat[gene.mat$'adj.P.Val' < 0.05,];
      gene.list$metadata <- list(gene=rownames(gene.mat),logFC=unname(gene.mat[,1]), adjP = unname(gene.mat$'adj.P.Val')); 
    } 
  }else{
    gene.mat <- dataSet$sig.mat;
    gene.list$metadata <- list(gene=rownames(gene.mat),logFC=unname(gene.mat$logFC), adjP = unname(gene.mat$'adj.P.Val'));
  }
  write.table(gene.mat, file="sig_genes.txt");
  gene.vec <- rownames(gene.mat);
  GeneAnotDB <- convertIdToEntrez(rownames(gene.mat), "entrez");
  protein.vec <- GeneAnotDB[,2];
  gene.mat <- data.matrix(gene.mat);
  rownames(gene.mat) <- protein.vec;
  na.inx <- is.na(protein.vec);
  prot.mat <- gene.mat[!na.inx,, drop=F];
  #write.table(cbind(Uniprot=rownames(prot.mat), Expression=prot.mat[,1]), file="seed_proteins.txt", row.names=F, quote=F);
  write.table(cbind(Emblprotein=rownames(prot.mat), Expression=prot.mat[,1]), file="seed_proteins.txt", row.names=F, quote=F);
  protein.vec <- prot.mat[,1];
  if(length(protein.vec) == 1){
    protein.vec <- as.matrix(protein.vec)
  }
  protein.list <- list();
  protein.list$metadata <- signif(protein.vec, 5);
  seed.expr <- prot.mat[,1];
  seed.df <- as.matrix(seed.expr);
  rownames(seed.df) <- names(seed.expr);
  seed.expr <- RemoveDuplicates(seed.df, "max", quiet=F); 
  seed.expr <<- seed.expr[,1];
  seed.genes <<- unique(gene.vec);
  protein.vec <- unique(names(protein.vec));
  
  list(
    gene.list = gene.list,
    protein.list = protein.list,
    protein.vec = protein.vec
  );
}


# read tab delimited file
# can have many classes, stored in meta.info (starts with #) 
# return a list (data.name, data.frame, meta.data)
.readTabData <- function(dataName) {
  if(length(grep('\\.zip$',dataName,perl=TRUE))>0){
    dataName <- unzip(dataName);
    if(length(dataName) > 1){
      # test if "__MACOSX" or ".DS_Store"
      osInx <- grep('MACOSX',dataName,perl=TRUE);
      if(length(osInx) > 0){
        dataName <- dataName[-osInx];
      }
      dsInx <- grep('DS_Store',dataName,perl=TRUE);
      if(length(dsInx) > 0){
        dataName <- dataName[-dsInx];
      }
      dat.inx <- grep(".[Tt][Xx][Tt]$", dataName);
      if(length(dat.inx) != 1){
        current.msg <<- "More than one text files (.txt) found in the zip file.";
        return(0);
      }
    }
  }
  
  msg <- NULL;
  # using the powerful fread function, 10 times faster, note: default return data.table, turn off
  dat1 <- .readDataTable(dataName);
  
  # look for #CLASS, could have more than 1 class labels, store in a list
  meta.info <- list();
  cls.inx <- grep("^#CLASS", dat1[,1]);
  if(length(cls.inx) > 0){ 
    for(i in 1:length(cls.inx)){
      inx <- cls.inx[i];
      cls.nm <- substring(dat1[inx, 1],2); # discard the first char #
      if(nchar(cls.nm) > 6){
        cls.nm <- substring(cls.nm, 7); # remove class
      }
      cls.lbls <- dat1[inx, -1];
      # test NA
      na.inx <- is.na(cls.lbls);
      cls.lbls[na.inx] <- "NA";
      cls.lbls <- ClearFactorStrings(cls.nm, cls.lbls);
      
      meta.info[[cls.nm]] <- cls.lbls;
    }
  }else{
    current.msg <<- "No metadata labels #CLASS found in your data!";
    return("F");
  }
  
  meta.info <- data.frame(meta.info);
  
  # now remove all comments in dat1
  # assign rownames after covert to matrix as data.frame does not allow duplicate names
  comments.inx <- grep("^#", dat1[,1]);
  dat1.nms <- dat1[-comments.inx,1];
  dat1<-dat1[-comments.inx,-1];
  dat1 <- data.matrix(dat1);
  rownames(dat1) <- dat1.nms;
  
  list(
    name= basename(dataName),
    data=dat1,
    meta.info=meta.info
  );
}


# note, try to use the fread, however, it has issues with 
# some windows 10 files "Line ending is \r\r\n. .... appears to add the extra \r in text mode on Windows"
# in such as, use the slower read.table method
.readDataTable <- function(fileName){
  if(length(grep('\\.zip$',fileName,perl=TRUE))>0){
    fileName <- unzip(fileName);
    if(length(fileName) > 1){
      # test if "__MACOSX" or ".DS_Store"
      osInx <- grep('MACOSX',fileName,perl=TRUE);
      if(length(osInx) > 0){
        fileName <- fileName[-osInx];
      }
      dsInx <- grep('DS_Store',fileName,perl=TRUE);
      if(length(dsInx) > 0){
        fileName <- fileName[-dsInx];
      }
      dat.inx <- grep(".[Tt][Xx][Tt]$", fileName);
      if(length(dat.inx) != 1){
        current.msg <<- "More than one text files (.txt) found in the zip file.";
        return(0);
      }
    }
  }
  dat <- try(data.table::fread(fileName, header=TRUE, check.names=FALSE, data.table=FALSE));
  if(class(dat) == "try-error"){
    # try to use "tr" to remove double return characters
    trFileName <- paste("tr -d \'\\r\' <", fileName);
    dat <- try(data.table::fread(trFileName, header=TRUE, check.names=FALSE, data.table=FALSE));
    if(class(dat) == "try-error"){
      print("Using slower file reader ...");
      formatStr <- substr(fileName, nchar(fileName)-2, nchar(fileName))
      if(formatStr == "txt"){
        dat <-try(read.table(fileName,header=TRUE,comment.char = "", check.names=F, as.is=T));
      }else{ # note, read.csv is more than read.table with sep=","
        dat <-try(read.csv(fileName,header=TRUE,comment.char = "", check.names=F, as.is=T));
      }  
    }
  }
  return(dat);
}

meanSdPlot <- function(x, ranks = TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
                       ylab = "sd", pch, plot = TRUE, bins = 50, ...) {
  
  stopifnot(is.logical(ranks), length(ranks) == 1, !is.na(ranks))
  
  n = nrow(x)
  if (n == 0L) {
    warning("In 'meanSdPlot': input matrix 'x' has 0 rows. There is nothing to be done.")
    return()
  }
  if (!missing(pch)) {
    warning("In 'meanSdPlot': 'pch' is ignored.")
  }
  
  px   = rowMeans(x, na.rm = TRUE)
  py   = sqrt(rowV(x, mean = px, na.rm = TRUE))
  rpx  = rank(px, na.last = FALSE, ties.method = "random")
  
  ## run median with centers at dm, 2*dm, 3*dm,... and width 2*dm
  dm        = 0.025
  midpoints = seq(dm, 1-dm, by = dm)
  within    = function(x, x1, x2) { (x >= x1) & (x <= x2) }
  mediwind  = function(mp) median(py[within(rpx/n, mp - 2*dm, mp + 2*dm)], na.rm = TRUE)
  rq.sds    = sapply(midpoints, mediwind)
  
  res = if(ranks) {
    list(rank = midpoints*n, sd = rq.sds, px = rpx, py = py)
  } else {
    list(quantile = quantile(px, probs = midpoints, na.rm = TRUE), sd = rq.sds, px = px, py = py)
  }
  
  fmt = function() function(x) format(round(x, 0), nsmall = 0L, scientific = FALSE)
  
  res$gg = ggplot(data.frame(px = res$px, py = res$py),
                  aes_string(x = "px", y = "py")) + xlab(xlab) + ylab(ylab) +
    geom_hex(bins = bins, ...) +
    scale_fill_gradient(name = "count", trans = "log", labels = fmt()) + 
    geom_line(aes_string(x = "x", y = "y"),
              data = data.frame(x = res[[1]], y = res$sd), color = "red")
  
  if (plot) print(res$gg)
  
  return(invisible(res))
}

rowV = function(x, mean, ...) {
  sqr     = function(x)  x*x
  n       = rowSums(!is.na(x))
  n[n<1]  = NA
  if(missing(mean))
    mean=rowMeans(x, ...)
  return(rowSums(sqr(x-mean), ...)/(n-1))
}

##################################################
## R script for NetworkAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# init resources for analysis
Init.Data<-function(path="../../"){
  selectedFactorInx <<- 1;
  fundbu <<- "kegg";
  chord_count <<-0;
  numOfLists <<- 1;
  rankOptGlobal <<- "pval";
  initmetaloading <<- TRUE;
  data.org <<- "hsa";
  keggpw.count <<- 0;
  pvalu <<- 0.05;
  dataSet <<- list(annotated=FALSE); 
  if(file.exists("/home/glassfish/sqlite/networkanalyst/")){
    sqlite.path <<- "/home/glassfish/sqlite/networkanalyst/";  #public server
    genesdb.path <<- "/home/glassfish/sqlite/"
    # disable parallel prcessing for DESeq2/edgeR
    library(BiocParallel);
    register(SerialParam());
  }else if(file.exists("/Users/xia/Dropbox/sqlite/networkanalyst/")){
    sqlite.path <<- "/Users/xia/Dropbox/sqlite/networkanalyst/"; #xia local
    genesdb.path <<- "/Users/xia/Dropbox/sqlite/"
  }else if(file.exists("/home/zzggyy/Downloads/netsqlite/")){
    sqlite.path <<-"/home/zzggyy/Downloads/netsqlite/"; #zgy local
    genesdb.path <<-"/home/zzggyy/Downloads/netsqlite/"
  }else if(file.exists("~/Downloads/sqlite")){
    sqlite.path <<- "~/Documents/Projects/sqlite/networkanalyst/"; #soufanom local
    genesdb.path <<- "~/Downloads/sqlite/"
  }
  
  lib.path <<- paste0(path, "data/");
  data.org <<- NULL; 
  module.count <<- 0; 
  msg.vec <<- vector(mode="character");
  current.msg <<- "";
  
  # preload some general package
  require('Cairo');
  CairoFonts("Arial:style=Regular","Arial:style=Bold","Arial:style=Italic","Helvetica","Symbol")
  require('igraph');
  print("called networkanalyst init!");
}

# genelist, onedata, metadata
# also set up or clear the other global objects
SetAnalType <- function(analType){
  anal.type <<- analType;
  mdata.all <<- list(); 
  meta.selected <<- TRUE;
  meta.upload <<- FALSE; # when upload merged data from meta-analysis b4
}

SetNetType <- function(netType){
  net.type <<- netType;
}

# When multiple genelists/datasets/results, record their name and save the data as .RDS file
# a) Current dataSet object
# Note, the memory will only contain one dataSet object. By default, the last one will be the current dataSet object;
# Users can switch this (from the interface) to specify which data is load into memory (dataSet object)
# b) Include for certain analysis
# For chord and heatmap analysis, users can do multiple selection (include)
# All datasets are selected by default (1 for selected, 0 for unselected)

# note, dataSet need to have "name" property
RegisterData <- function(dataSet){
  dataName <- dataSet$name;
  saveRDS(dataSet, file=dataName);
  if(!is.null(dataSet$data.raw)){ # save memory for meta-analysis mode
    dataSet$data.raw <- NULL;
  }
  dataSet <<- dataSet;
  if(anal.type == "metadata"){
    mdata.all[[dataName]] <<- 1;
  }else{
    mdata.all <<- lapply(mdata.all, function(x){ x <- 0;});
    mdata.all[[dataName]] <<- 1;
  }
  return(1);
}

# only for switching single expression data results
SetCurrentData <- function(nm){
  if(dataSet$name != nm){
    dataSet <<- readRDS(nm);
  }
  return(1);
}

# remove data object, the current dataSet will be the last one by default 
RemoveData <- function(dataName){
  if(!is.null(mdata.all[[dataName]])){
    mdata.all[[dataName]] <<- NULL;
  }
}

# users can select one or more data for analysis
# note, we use 1 to indicate this is selected
# and by default is all selected. 
SelectData <- function(){
  if(!exists('nm.vec')){
    current.msg <<-"No dataset is selected for analysis!";
    print(current.msg);
    return(0);
  }
  
  all.nms <- names(mdata.all);
  for(nm in all.nms){
    if(nm %in% nm.vec){
      mdata.all[[nm]] <<- 1;
    }else{
      mdata.all[[nm]] <<- 0;
    }
  }
  if(anal.type == "metadata"){
    if("meta_dat" %in% nm.vec){
      meta.selected <<- TRUE;
    }else{
      meta.selected <<- FALSE;
    }
  }
  rm('nm.vec', envir = .GlobalEnv);
  return(1);
}

GetAllDataNames <- function(){
  names(mdata.all);
}

SetOrganism <- function(org){
  data.org <<- org;
  init.lib <<- "kegg"
}

SetSelectedFactorInx <- function(inx){
  selectedFactorInx <<- inx;
}

SetSelNetDataset <- function(type){
  selectedNetDataset <<- type;
}

SetRankingMetric <- function(opt){
  rankOptGlobal <<- opt;
}


SetListNms <- function(){
  newDat <- list();
  tot.count <- 0;
  listSizes <- list();
  
  # convert to entrez
  if(anal.type == "metadata"){
    inmex.meta <- readRDS("inmex_meta.rds");
    en.ids <- rownames(inmex.meta$data);
    nm = "meta_data"
  }else{
    en.ids <- rownames(dataSet$resTable)
    nm = "dataSet"
  }
  names(en.ids) <- doEntrez2SymbolMapping(en.ids)
  
  listSizes[[1]] <- list(
    name = nm,
    label = nm,
    size = length(en.ids)
  );
  
  list.genes <<- en.ids;
  listSizes <<- listSizes;
}

SetInitLib <- function(org){
  init.lib <<- "kegg"
}

GetDataListNames <- function(){
  return(names(mdata.all));
}

##################################################
## R script for NetworkAnalyst
## Description: functions only for list data analysis
##
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

SaveListHeatmapJSON <- function(fileName){
  if(numOfLists>1){
    json.res <- PrepareMultiListHeatmapJSON();
  }else{
    json.res <- PrepareListHeatmapJSON();
  }
  require(RJSONIO);
  json.mat <- toJSON(json.res, .na='null');
  sink(fileName);
  cat(json.mat);
  sink();
  current.msg <<- "Data is now ready for heatmap visualization!";
  return(1);
}

PrepareListHeatmapJSON <- function(){
  sig.ids <- rownames(dataSet$prot.mat);
  gene.symbols=doEntrez2SymbolMapping(sig.ids)
  stat.pvals <- dataSet$prot.mat[,1]
  
  expval <- 0
  expval <- sum(dataSet$prot.mat)
  
  # scale each gene 
  #data.stat <- readRDS("data.stat");
  #hit.inz = sig.ids %in% rownames(data.stat);
  #sig.ids = sig.ids[hit.inz];
  dat <- dataSet$prot.mat
  
  # now pearson and euclidean will be the same after scaleing
  dat.dist <- dist(dat); 
  
  orig.smpl.nms <- colnames(dat);
  orig.gene.nms <- rownames(dat);
  
  # prepare meta info    
  # 1) convert meta.data info numbers
  # 2) match number to string (factor level)
  
  grps <- "datalist1"
  cls <- "datalist1"
  
  # convert back to numeric 
  
  # for each gene/row, first normalize and then tranform real values to 30 breaks
  if(expval !=0){
    dat_pos = as.matrix(dat[sign(dat[,1]) == 1,])
    dat_neg = as.matrix(dat[sign(dat[,1]) == -1,])
    if(nrow(dat_pos) == 0){
      res <- apply(unname(dat), 2, function(x){
        y =log(abs(x)) + 0.000001
        16-as.numeric(cut(y, breaks=15))
      });
    }else if(nrow(dat_neg) == 0){
      res <- apply(unname(dat), 2, function(x){
        y =log(x) + 0.000001
        15+as.numeric(cut(y, breaks=15))
      });
    }else{
      res_pos <- apply(unname(dat_pos), 2, function(x){
        y =log(x) + 0.000001
        as.numeric(cut(y, breaks=15))+15
      });
      res_neg <- apply(unname(dat_neg), 2, function(x){
        y =log(abs(x)) + 0.000001
        16 - as.numeric(cut(y, breaks=15))
      });
      res = rbind(res_pos, res_neg);
    }
  }else{
    zero.inx = dataSet$prot.mat == 0
    res = dataSet$prot.mat;
    res[zero.inx] = 32
  }
  
  res_list <- list()
  for(i in 1:length(res)){
    res_list[[i]] <- list(res[i])
  }
  
  # note, use {} will lose order; use [[],[]] to retain the order
  
  nmeta = list(100)
  nmeta.anot = list()
  
  nmeta.anot["datalist1"] = nmeta[1]
  
  nmeta = list(nmeta)
  names(nmeta) = "datalists"
  
  json.res <- list(
    data.type = "singlelist", 
    gene.id = gene.symbols,
    gene.entrez = sig.ids,
    gene.name = gene.symbols,
    gene.cluster = 1,
    sample.cluster = 1,
    sample.names = list("datalist1"),
    meta = nmeta,
    meta.anot = nmeta.anot,
    data = res_list,
    expval = expval
  );
  rownames(dat) = gene.symbols
  write.csv(dat,"heatmap_matrix.csv", row.names=TRUE)
  return(json.res);
}

PrepareMultiListHeatmapJSON <- function(){
  sel.nms <- names(mdata.all)
  expval<-0;
  for(i in 1:length(sel.nms)){
    dataNm <- sel.nms[i];
    dataSet <- readRDS(dataNm);
    len <- nrow(dataSet$prot.mat)
    if(i == 1){
      expval <- sum(dataSet$prot.mat)
      gene_list <-rownames(dataSet$prot.mat)
    }else{
      gene_list <-c(gene_list, rownames(dataSet$prot.mat))
      expval <- expval + sum(dataSet$prot.mat)
    }
  }
  
  gene_list <- unique(gene_list)
  allmat = matrix(NA, nrow=length(gene_list), ncol=length(sel.nms))
  rownames(allmat) = gene_list
  
  for(i in 1:length(sel.nms)){
    dataName <- sel.nms[i];
    dataSet <- readRDS(dataName);
    cols <- colnames(allmat)[colnames(allmat) %in% dataName]
    if(expval ==0){
      rows <- which(rownames(allmat) %in% rownames(dataSet$prot.mat))
      inx <-match(rownames(allmat) ,rownames(dataSet$prot.mat))
      allmat[, i] <- as.vector(dataSet$prot.mat)[inx]
    }else{
      rows <- which(rownames(allmat) %in% rownames(dataSet$prot.mat))
      inx <-match(rownames(allmat) ,rownames(dataSet$prot.mat))
      allmat[, i] <- as.vector(dataSet$prot.mat)[inx]
    }
  } 
  colnames(allmat) = sel.nms 
  inx <- apply(allmat, 1, function(x){sum(is.na(x))});  
  ord.inx <- order(inx)
  allmat = allmat[ord.inx,]
  gene.symbols = doEntrez2SymbolMapping(rownames(allmat))
  
  na.inx = is.na(allmat)
  zero.inx = allmat == 0
  
  allmatb = allmat
  
  allmatb[na.inx]=0
  allmatb[zero.inx]=1
  rownames(allmatb) = gene.symbols
  write.csv(allmatb,"heatmap.csv", row.names=TRUE)  
  
  if(expval != 0){
    pos.inx = allmat>0 & !na.inx
    neg.inx = allmat<0 & !na.inx
    allmat[neg.inx] = 16 - as.numeric(cut(log(abs(allmat[neg.inx])) , breaks=15))
    allmat[pos.inx] = 15 + as.numeric(cut(log(allmat[pos.inx]) , breaks=15))
    allmat[zero.inx] = 32
  }else{
    zer.inx = allmat == 0 & !na.inx
    nb <- apply(allmat, 1, function(x){sum(!is.na(x))});
    for(i in 1:nrow(allmat)){
      row = allmat[i,]
      inx = row == 0
      allmat[i, inx] = nb[i]
    }
    allmat[zer.inx] <- try(15 + as.numeric(cut(allmat[zer.inx] , breaks=15)));
    if(class(allmat[zer.inx]) == "try-error") {
      allmat[zer.inx] = 32
    }else{
      15 + as.numeric(cut(allmat[zer.inx] , breaks=15))
    }
  }
  allmat[na.inx] = 31;
  res_list <- list();
  for(i in 1:nrow(allmat)){
    res_list[i] <- list(unname(allmat[i,]))
  }
  
  nmeta = as.numeric(as.factor(colnames(allmat))) + 99
  nmeta.anot = list()
  
  for(i in 1:length(unique(nmeta))){
    nmeta.anot[[colnames(allmat)[i]]] = nmeta[i]
  }
  nmeta = list(nmeta)
  names(nmeta) = "datalists"
  
  json.res <- list(
    data.type = "mutlilist",
    gene.id = gene.symbols,
    gene.entrez = rownames(allmat),
    gene.name = rownames(allmat),
    gene.cluster = 1,
    sample.cluster = 1,
    sample.names = colnames(allmat),
    meta = nmeta,
    meta.anot = nmeta.anot,
    data.lbl = "NA",
    data = res_list,
    expval = expval
  );
  
  return(json.res);
}

##################################################
## R script for NetworkAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# read individual data from user, data.type: array
ReadIndData <- function(dataName, data.type){
  
  current.msg <<- "";
  meta.upload <<- FALSE; # upload data to resume
  dataSet <- .readTabData(dataName);
  
  # now need to remove low quality samples and genes
  data <- dataSet$data;
  meta.info <- dataSet$meta.info;
  
  smpl.num <- ncol(data);
  gene.num <- nrow(data);
  
  # remove smpls/exp with over half missing value
  good.inx<-apply(is.na(data), 2, sum)/nrow(data)<0.6;
  smpl.msg <- "";
  if(sum(!good.inx)>0){
    smpl.msg <- paste(sum(!good.inx), "low quality samples(>60% missing) removed.");
    print(smpl.msg);
    
    data <- data[,good.inx];
    if(ncol(data)/smpl.num < 0.5){
      current.msg <<- paste(smpl.msg, "Low quality data rejected!");;
      return("F");
    }
    
    # update meta information
    meta.info <- meta.info[good.inx, , drop=F];
  }
  
  if(ncol(data) < 4){
    current.msg <<- paste(smpl.msg, "The sample # (", ncol(data), ") is too small.");
    return("F");
  }
  
  # genes with 75% NA will be removed
  gd.inx<-apply(is.na(data), 1, sum)/ncol(data)<0.75;
  feat.msg <- "";
  if(sum(!gd.inx) > 0){
    data <- data[gd.inx,];
    feat.msg <- paste(sum(!gd.inx), "low quality genes (>75% missing) removed");
    if(nrow(data)/gene.num < 0.25){
      current.msg <<- paste(feat.msg, "Low quality data rejected.");
      return("F");
    }
    print(feat.msg);
  }
  
  if(nrow(data) < 10){ 
    current.msg <<- paste(feat.msg, "The gene# (", nrow(data), ") is too small (<10).");
    return("F");
  }
  
  # make an copy, only for testing different normalization
  dataSet$data.raw <- data;
  dataSet$data <- data;
  dataSet$type <- data.type;
  dataSet$meta.info <- meta.info;
  dataName <- dataSet$name;
  res <- RegisterData(dataSet);
  if(res == 1){
    return(dataName);
  }else{
    current.msg <<- paste("Cannot add data: ", dataName, ". ", current.msg, sep="");
    return("F");
  }
}

# read in the data and perform
# gene ID mapping using built in libraries
# matchMin is minimal matched probe (%)
# return the total matched gene number
ProcessIndData<-function(dataName, featureType, matchMin=0.5){
  
  dataSet <- readRDS(dataName);
  dataSet$id.type <- featureType;
  if(data.org != 'NA' & featureType != 'NA'){
    feature.vec <- rownames(dataSet$data.raw);
    minLvl <- length(feature.vec)*matchMin;
    if(featureType %in% c("entrez", "symbol", "refseq", "genbank", "emblgene", "embltranscript", "orfid", "wormbase")){
      entrez.id <- doGeneIDMapping(feature.vec, featureType);
    }else{
      entrez.id <- doProbeMapping(feature.vec, featureType);
    }   
    
    hit.inx <- which(!is.na(entrez.id));
    matched.len <- length(hit.inx);
    if(matched.len > minLvl){
      res <- "Success!";
      data.orig <- dataSet$data.raw[hit.inx,];
      matched.entrez <- entrez.id[hit.inx];
      rownames(data.orig) <- matched.entrez;
      
      # now, deal with duplicated entrez id
      # first, average duplicate rows
      ave.data <- apply(data.orig, 2, myave, matched.entrez); 
      # then removed duplicated entries
      dup.inx <- duplicated(matched.entrez);
      int.mat <- ave.data[!dup.inx,];
      # update
      dataSet$data.orig <- int.mat;
      dataSet$id.type <- "entrez";
    }
    current.msg <<- paste("ID Conversion: ", "Total [", length(entrez.id), 
                          "] Matched [", matched.len, "] Unmatched [", sum(is.na(entrez.id)),"]", collapse="\n"); 
  }else{ # no conversion will be performed
    dataSet$data.orig <- dataSet$data.raw;
    matched.len <- 9; # dummies
    minLvl <- 1;
  }
  RegisterData(dataSet);
  if(matched.len < minLvl){
    return(0);
  }
  return(matched.len);
}

# overwrite ave, => na.rm=T
myave <- function (x, ...) {
  n <- length(list(...))
  if (n) {
    g <- interaction(...)
    split(x, g) <- lapply(split(x, g), mean, na.rm=T)
  }
  else x[] <- FUN(x, na.rm=T)
  x
}

# read the uploaded data into memory
# return the meta-data information (multiple groups)
ReadDataForMetaInfo<-function(dataName){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  return(colnames(dataSet$meta.info));
}

# here should first try to load the original data
# the data in the memory could be changed
GetGroupNames <- function(dataName){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  levels(dataSet$cls);
}

GetDataDims <- function(dataName){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  dm <- dim(dataSet$data);
  naNum <- sum(is.na(dataSet$data));
  return(c(dm, naNum));
} 

# obtain sample names and their class labels
GetSampleInfo <- function(dataName, clsLbl){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  grpInfo <- dataSet$meta.info[[clsLbl]];
  grpLbls <- paste(levels(grpInfo), collapse="\n");
  smplInfo <- paste(Sample = colnames(dataSet$data.orig), "\t", Class=grpInfo, collapse="\n");
  return(c(grpLbls, smplInfo));
}

# given dataSet Name, sample name, and class name, do update
# note, for multiple #class, this set which one to use in the subsequent steps
# last one wins
UpdateSampleInfo<-function(dataName, clsLbl){
  
  print("updating sample info .... ");
  
  if(!exists("class.vec")){
    print("Could not find class label list!");
    return(0);
  }
  
  if(!exists("smpl.vec")){
    print("Could not find sample name list!");
    return(0);
  }
  
  if(length(class.vec) < 2){
    current.msg <<- "Add least two groups required!";
    return(0);
  }
  
  if(sum(class.vec != 'NA') < 2){
    current.msg <<- "Cannot be less than 2 groups";
    return(0);
  }
  
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  org.lvl.len <- length(levels(dataSet$meta.info[[clsLbl]]));
  if(org.lvl.len < length(class.vec)){
    current.msg <<- "You can not add new groups";
    return(0);
  }else if(org.lvl.len > length(class.vec)){
    current.msg <<- "To exclude a group, replace it with NA.";
    return(0);
  }
  
  # first update the meta info
  cls <- dataSet$meta.info[[clsLbl]];
  levels(cls) <- class.vec;
  
  data <- dataSet$data.orig;
  meta.info <- dataSet$meta.info;
  
  if(any(levels(cls) == 'NA')){
    rt.inx <- cls != 'NA';
    data <- data[,rt.inx];
    
    # also update the whole meta-info
    meta.info <- meta.info[rt.inx,,drop=FALSE];
    cls <- cls[rt.inx];
  }
  
  # need to re-construct the class, so that the level order  
  # are always alphabetic
  meta.info[[clsLbl]] <- factor(as.character(cls));
  
  # note, sample names could be removed (together with cls) as the whole row
  hit.inx <- colnames(data)%in%smpl.vec;
  dataSet$data.orig <- data[,hit.inx];
  
  # make sure the factor levels also dropped
  for(i in 1:length(meta.info)){
    meta.info[[i]] <- factor(meta.info[[i]][hit.inx]);
  }
  
  dataSet$meta.info <- meta.info;
  dataSet$cls <-  dataSet$meta.info[[clsLbl]];
  RegisterData(dataSet);
  gc();
  return(1);
}


# note, here also update data type array/count
PerformIndNormalization <- function(dataName, norm.opt, auto.opt, dataType){
  
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  msg <- NULL;
  data <- dataSet$data.orig;
  data <- PerformDataNormalization(data, norm.opt);
  if(length(data)==1 && data == 0){
    return(0);
  }
  msg <- paste(norm.msg, msg);
  
  if(auto.opt==1){
    row.nms <- rownames(data);
    col.nms <- colnames(data);
    data<-apply(data, 2, AutoNorm);
    msg <- paste(msg, "Autoscaling performed.", collapse=" ");
    rownames(data) <- row.nms;
    colnames(data) <- col.nms;
  }
  
  dataSet$data <- data;
  dataSet$type <- dataType;
  RegisterData(dataSet);
  current.msg <<- msg;
  return(1);
}

# normalize to zero mean and unit variance
AutoNorm<-function(x){
  (x - mean(x))/sd(x, na.rm=T);
}

######################################
## methods for merged expression data
#######################################

GlobalCutOff = list(
  logFC = 0,
  BHth = 0.05
)

# read meta-dataset previously processed
ReadMergedExpressTable <- function(dataName){
  current.msg <<- "";
  meta.upload <<- TRUE;
  dataSet <- .readTabData(dataName);
  common.matrix <- dataSet$data;
  meta.nms <- tolower(names(dataSet$meta.info));
  
  cls.inx <- grep("condition", meta.nms);
  if(length(cls.inx) == 0){
    current.msg <<- "No condition label found (#CLASS.condition)";
    return("F");
  }else{
    cls.inx <- cls.inx[1];
    cls.lbl <- dataSet$meta.info[[cls.inx]];
  }
  
  data.inx <- grep("dataset", meta.nms);
  if(length(data.inx) == 0){
    current.msg <<- "No dataset label found (#CLASS.dataset)";
    return("F");
  }else{
    data.inx <- data.inx[1];
    data.lbl <- dataSet$meta.info[[data.inx]];
    data.nms <- unique(as.character(data.lbl));
    
    # now create the mdata.all object
    mdata.all <- vector(mode="list", length=length(data.nms));
    names(mdata.all) <- data.nms;
    mdata.all <<- lapply(mdata.all, function(x){x=1});
  }
  
  if(length(grep("entrez.hsa", meta.nms)) > 0){
    data.org <<- "hsa"
    id.type <<- "entrez";
    shared.nms <- rownames(common.matrix);
    symbols <- doEntrez2SymbolMapping(shared.nms);
    names(symbols) <- shared.nms;
  }else if(length(grep("entrez.mmu", meta.nms)) > 0){
    data.org <<- "mmu"
    id.type <<- "entrez";
    shared.nms <- rownames(common.matrix);
    symbols <- doEntrez2SymbolMapping(shared.nms);
    names(symbols) <- shared.nms;
  }else{
    symbols <- NULL;
    inmex.org <<- "NA"
    id.type <<- 'NA';
  }
  data.org <<- unlist(strsplit( meta.nms[2], "[.]"))[3]
  inmex.meta.orig <- list(data=common.matrix,
                          id.type = id.type,
                          gene.symbls = symbols,
                          cls.lbl=factor(cls.lbl),
                          data.lbl=data.lbl);
  saveRDS(inmex.meta.orig, "inmex.meta.orig.rds");
  if(length(levels(as.factor(data.lbl))) == 1){
    return(2);
  }else{
    return(1);
  }    
}

GetMetaMetaInfo <- function(dataName){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  return(colnames(dataSet$meta.info));
}

SetSelectedMetaMetaInfo <- function(dataName, meta0, meta1, block1){
  if(meta0 == "NA"){
    return(0);
  }else{
    if(dataSet$name != dataName){
      dataSet <- readRDS(dataName);
    }
    cls <- dataSet$meta.info[, meta0];
    block <- NULL;
    if(meta1 != "NA"){
      if(block1){
        block <- dataSet$meta.info[, meta1]
      }else{ # two factor
        cls <- interaction(dataSet$meta.info[, c(meta0, meta1)], sep = ".", lex.order = TRUE);
      }
    }
    dataSet$cls <- cls; # record main cls;
    dataSet$block <- block;
    RegisterData(dataSet);
    gc();
    return(levels(cls));
  }
}


# for gene-level meta-analysis
# function to set up results combining individual data analysis
# as well as to prepare for GO analysis
# no return, as set global 

SetupMetaStats <- function(BHth){
  
  GlobalCutOff$BHth <<- BHth;
  #all common genes
  inmex.meta <- readRDS("inmex_meta.rds");
  gene.ids <- rownames(inmex.meta$data);
  # meta.sig genes
  metade.genes <- rownames(meta.mat);
  
  # setup individual sig genes & stats
  # that overlap with meta.sig
  inmex.de <- list();
  
  pval.mat <- fc.mat <- matrix(nrow=nrow(meta.mat), ncol=sum(mdata.all==1));
  
  for(i in 1:length(inmex.ind)){
    de.res <- inmex.ind[[i]];
    
    hit.inx <- de.res[,2] <= BHth;
    hit.inx <- which(hit.inx); # need to get around NA
    inmex.de[[i]] <- rownames(de.res)[hit.inx];
    
    # only choose the genes that are also meta sig genes from in
    # individual analysis for display
    de.res <- de.res[metade.genes,];
    
    fc.mat[,i] <- de.res[,1];
    pval.mat[,i] <- de.res[,2];
  }
  names(inmex.de) <- names(inmex.ind);
  
  # calculate gain/loss
  deindst <- unique(unlist(inmex.de));
  gains=metade.genes[which(!(metade.genes %in% deindst))];
  losses=deindst[which(!(deindst %in% metade.genes))];
  all.de <- cbind(gene.ids %in% metade.genes, gene.ids %in% deindst);
  colnames(all.de) <- c("Meta-DE", "Individual-DE");
  vennC <- getVennCounts(all.de);
  if(inmex.meta$id.type == "entrez"){ 
    names(metade.genes) <- inmex.meta$gene.symbls[metade.genes];
    names(gains) <- inmex.meta$gene.symbls[gains];
    names(losses) <- inmex.meta$gene.symbls[losses];
  }
  
  # de genes from individual 
  de.len <- sapply(inmex.de, length);
  stat <- c(length(metade.genes), de.len);
  names(stat) <- c("Meta", substr(names(inmex.de), 0, nchar(names(inmex.de))-4));
  meta.stat <- list(
    stat = stat,
    de = metade.genes,
    idd = gains,
    loss = losses,
    venn = vennC
  );
  
  fc.mat <<- fc.mat;
  pval.mat <<- pval.mat;
  inmex.de <<- inmex.de;
  meta.stat <<- meta.stat;
  
  # save the result
  if(inmex.meta$id.type == "entrez"){ # row name gene symbols
    metade.nms <- inmex.meta$gene.symbls[metade.genes];
    res <- cbind(EntrezID=metade.genes, Name=metade.nms, meta.mat);
  }else{
    res <- cbind(ID=metade.genes, meta.mat);
  }
  write.csv(res, file=paste("meta_sig_genes_", inmex.method, ".csv", sep=""), row.names=F);
}

# perform differential analysis
# default: all pair-wise comparison (A-B) + (B-C) + (A-C)
# custom: only compare two groups (A-C)
# time: compare consecutive groups (B-A) + (C-B)
# reference: all others against common reference (A-C) + (B-C)
# nested: (A-B)+(C-D) 
PerformIndLimmaDE<-function(dataName, anal.type, par1=NULL, par2=NULL, nested.opt = "intonly"){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  require(limma);
  current.msg <<- "";
  cls <- dSet$cls; 
  
  cls.lvls <- levels(cls);
  data <- dSet$data.norm;
  meta.info <- dSet$meta.info;
  
  design <- model.matrix(~ -1 + cls) # no intercept
  colnames(design) <- cls.lvls;
  
  if(is.null(dSet$block)){
    fit = lmFit(data, design);
  }else{
    # limma user guide chapter 8 P49 
    corfit<-duplicateCorrelation(data,design,block=dSet$block);
    fit = lmFit(data, design, block=dSet$block, correlation=corfit$consensus);
  }
  # sanity check
  if(!is.fullrank(design)){
    current.msg <<- paste("This metadata combination is not full rank! Please use other combination."); 
    return(0);
  }
  df.residual <- fit$df.residual;
  if (all(df.residual == 0)){
    current.msg <<- paste("There is not enough replicates in each group (no residual degrees of freedom)!"); 
    return(0);
  }
  
  myargs <- list();
  grp.nms <- cls.lvls;
  
  if(anal.type == 'default'){ # all pair-wise
    inx = 0;
    for(m in 1:(length(grp.nms)-1)){
      for(n in (m+1):length(grp.nms)){
        inx <- inx + 1;
        myargs[[inx]] <- paste(grp.nms[m], "-", grp.nms[n], sep="")
      }
    }
  }else if(anal.type == 'time'){
    for(i in 2:length(grp.nms)){
      myargs[[i-1]] <- paste(grp.nms[i], "-", grp.nms[i-1], sep="")
    }
  }else if(anal.type == 'custom'){
    grp.nms <- strsplit(par1, " vs. ")[[1]];
    myargs[[1]] <- paste(grp.nms, collapse="-");
  }else if(anal.type == 'reference'){
    ref <- par1;
    cntr.cls <- grp.nms[grp.nms!= ref];
    myargs <- as.list(paste(cntr.cls, "-", ref, sep=""));
  }else if(anal.type == 'nested'){
    grp.nms1 <- strsplit(par1, " vs. ")[[1]];
    grp.nms2 <- strsplit(par2, " vs. ")[[1]];
    if(all(grp.nms1 == grp.nms2)){
      current.msg <<- paste("The two nested groups are the same. Please choose two different groups."); 
      return(0);
    }
    grp.nms <- unique(c(grp.nms1, grp.nms2));
    if(nested.opt == "intonly"){
      myargs[[1]] <- paste("(", paste(grp.nms1, collapse="-"), ")-(", paste(grp.nms2, collapse="-"), ")", sep=""); 
    }else{
      myargs[[1]] <- paste(grp.nms1, collapse="-");
      myargs[[2]] <- paste(grp.nms2, collapse="-"); 
      myargs[[3]] <- paste("(", paste(grp.nms1, collapse="-"), ")-(", paste(grp.nms2, collapse="-"), ")", sep=""); 
    }
  }else{ # 
    print(paste('Not supported: ', anal.type));
  }
  
  myargs[["levels"]] <- design;
  contrast.matrix <- do.call(makeContrasts, myargs);
  fit2 <- contrasts.fit(fit, contrast.matrix);
  fit2 <- eBayes(fit2);
  
  resTable <- topTable(fit2, number=Inf, adjust.method="fdr");
  dSet$sig.orig <- resTable; # record orignal data for update
  anot.id <- rownames(resTable);
  gene.anot <- doEntrezIDAnot(anot.id);
  current.msg <<- current.msg; 
  
  # may need to update data and class
  hit.inx <- cls %in% grp.nms;
  if(sum(hit.inx) < length(hit.inx)){
    current.msg <<- paste(current.msg, "Only groups selected for comparisons: ", paste(grp.nms, collapse=", "), "are included.");
    cls <- factor(cls[hit.inx]);
    cls.lvls <- levels(cls);
    data <- data[,hit.inx];
    meta.info <- dSet$meta.info[hit.inx,];
  }
  
  dSet$data.stat <- data;
  dSet$cls.stat <- cls;
  dSet$meta.stat <- meta.info;
  dSet$sig.genes.anot <- gene.anot;
  RegisterData(dataSet);
  return (1);
}

GetProjectOrg <- function(){
  data.org;
}

GetProjectIDType <- function(){
  id.type;
} 

#  VENN DIAGRAM COUNTS AND PLOTS
getVennCounts <- function(x,include="both") {
  x <- as.matrix(x)
  include <- match.arg(include,c("both","up","down"))
  x <- sign(switch(include,
                   both = abs(x),
                   up = x > 0,
                   down = x < 0
  ))
  nprobes <- nrow(x)
  ncontrasts <- ncol(x)
  names <- colnames(x)
  if(is.null(names)) names <- paste("Group",1:ncontrasts)
  noutcomes <- 2^ncontrasts
  outcomes <- matrix(0,noutcomes,ncontrasts)
  colnames(outcomes) <- names
  for (j in 1:ncontrasts)
    outcomes[,j] <- rep(0:1,times=2^(j-1),each=2^(ncontrasts-j))
  xlist <- list()
  for (i in 1:ncontrasts) xlist[[i]] <- factor(x[,ncontrasts-i+1],levels=c(0,1))
  counts <- as.vector(table(xlist))
  structure(cbind(outcomes,Counts=counts),class="VennCounts")
}

#compute Cochran’s Q to help FEM/REM
# plot Q-Q plot for estimation
PlotCochranQ <- function(imgNm){
  if(!performedDE){
    PerformMetaDeAnal();
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  nbstudies <- sum(sel.inx);
  ES<-array(dim=c(nrow(inmex.meta$data),4,nbstudies));
  cls.lvls <- levels(as.factor(inmex.meta$cls.lbl));
  
  for (i in 1:nbstudies){
    data.nm <- sel.nms[i];
    dataSet <- readRDS(data.nm);
    
    fit.obj.nm <- paste(data.nm, "fit.obj", sep=".");
    fit2i <- readRDS(fit.obj.nm);
    #fit2i <- dataSet$fit.obj;
    
    n1i=length(which(dataSet$cls==cls.lvls[1]));
    n2i=length(which(dataSet$cls==cls.lvls[2]));
    ES[,,i]=effectsize(fit2i$t,((n1i*n2i)/(n1i+n2i)),(fit2i$df.prior+fit2i$df.residual));
  }
  
  Qvals <- f.Q.NA(ES[,3,],ES[,4,]);
  
  Cairo(file=imgNm, width=400, height=400, type="png", bg="white");
  # histgram
  # hist(Qvals, breaks = 50, col = "red");
  # QQ plot
  chisqq <- qchisq(seq(0, 0.9999, 0.001), df = nbstudies-1)
  tmp <- quantile(Qvals, seq(0, 0.9999, 0.001))
  qqplot(chisqq, tmp, ylab = "Quantiles of Sample", pch = "*", 
         xlab = "Quantiles of Chi > square", main = "QQ Plot")
  lines(chisqq, chisqq, lty = "dotted", col = "red")
  
  dev.off(); 
  
  Qmean <- round(mean(Qvals),5);
  return (Qmean);
}

PlotDataProfile<-function(dataName, boxplotName, pcaName){
  dataSet <- readRDS(dataName);
  qc.boxplot2(dataSet$data, boxplotName);
  qc.pcaplot2(dataSet$data, pcaName);
}


qc.boxplot2 <- function(dat, imgNm){
  require('lattice');
  imgNm = paste(imgNm, "dpi", "72", ".png", sep="");
  subgene=10000;
  if (nrow(dat)>subgene) {
    set.seed(28051968);
    sg  = sample(nrow(dat), subgene)
    Mss = dat[sg,,drop=FALSE]
  } else {
    Mss = dat
  }
  
  subsmpl=100;
  if (ncol(Mss)>subsmpl) {
    set.seed(28051968);
    ss  = sample(ncol(Mss), subsmpl)
    Mss = Mss[,ss,drop=FALSE]
  } else {
    Mss = Mss
  }
  
  sample_id = rep(seq_len(ncol(Mss)), each = nrow(Mss));
  values  = as.numeric(Mss)
  formula = sample_id ~ values
  
  box = bwplot(formula, groups = sample_id, layout = c(1,1), as.table = TRUE,
               strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
               horizontal = TRUE,
               pch = "|",  col = "black", do.out = FALSE, box.ratio = 2,
               xlab = "", ylab = "Samples",
               fill = "#1c61b6AA",
               panel = panel.superpose,
               scales = list(x=list(relation="free"), y=list(axs="i")),
               ylim = c(ncol(Mss)+0.7,0.3),
               prepanel = function(x, y) {
                 list(xlim = quantile(x, probs = c(0.01, 0.99), na.rm=TRUE))
               },
               panel.groups = function(x, y, ...) {
                 panel.bwplot(x, y, ...)
               })
  
  Cairo(file=imgNm, width=460, height=420, type="png", bg="white");
  print(box);
  dev.off();
}

qc.pcaplot2 <- function(x, imgNm){
  imgNm = paste(imgNm, "dpi", "72", ".png", sep="");
  require('lattice');
  pca <- prcomp(t(na.omit(x)));
  names <- colnames(x);
  pca.res <- as.data.frame(pca$x);
  # increase xlim ylim for text label
  xlim <- GetExtendRange(pca.res$PC1);
  ylim <- GetExtendRange(pca.res$PC2);
  pcafig = xyplot(PC2~PC1, data=pca.res, pch=19, cex=1,aspect = "iso", xlim = xlim, ylim=ylim,
                  panel=function(x, y, ...) {
                    panel.xyplot(x, y, ...);
                    ltext(x=x, y=y, labels=names, pos=1, offset=1, cex=0.8, col="magenta")
                  })
  
  Cairo(file=imgNm, width=480, height=480, type="png", bg="white");
  print(pcafig);
  dev.off();
}

PlotMetaPCA <- function(imgNm, dpi, format,factor){
  inmex.meta <- readRDS("inmex_meta.rds");
  x=inmex.meta[["data"]]
  dpi = as.numeric(dpi);
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  library('lattice');
  library('ggplot2');
  pca <- prcomp(t(na.omit(x)));
  imp.pca<-summary(pca)$importance;
  xlabel = paste0("PC1"," (", 100*round(imp.pca[2,][1], 3), "%)")
  ylabel = paste0("PC2"," (", 100*round(imp.pca[2,][2], 3), "%)")
  names <- colnames(x);
  pca.res <- as.data.frame(pca$x);
  # increase xlim ylim for text label
  xlim <- GetExtendRange(pca.res$PC1);
  ylim <- GetExtendRange(pca.res$PC2);
  if(factor != "NA"){
    #Factor = as.vector(dataSet$meta.info[,factor])
  }else{
    #Factor = dataSet$meta.info[,1];
  }
  Conditions = factor(inmex.meta$cls.lbl)
  Datasets = factor(inmex.meta$data.lbl)
  pcafig = ggplot(pca.res, aes(x=PC1, y=PC2,  color=Conditions ,shape=Datasets)) +
    geom_point(size=4, alpha=0.5) + xlim(xlim)+ ylim(ylim) + xlab(xlabel) + ylab(ylabel);
  
  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", unit="in", dpi=dpi);
  print(pcafig);
  dev.off();
  
}


GetMetaSummary<- function(){
  inmex.meta <- readRDS("inmex_meta.rds");
  sel.nms <- unique(inmex.meta$data.lbl)
  sel.nms <- paste(sel.nms, collapse="; ")
  cls.lbls <- unique(inmex.meta$cls.lbl)
  cls.lbls <- paste(cls.lbls, collapse="; ")
  return(c(length(colnames(inmex.meta$data)),nrow(inmex.meta$data), sel.nms, cls.lbls))
}

GetMetaDatasets<- function(){
  sel.nms <- names(mdata.all)[mdata.all==1];
  return(sel.nms);
}

SetSelMetaData<- function(selNm){
  selDataNm <<- selNm
}

# retrun the json obj
SaveMetaClusterJSON <- function(fileName, clustOpt, opt){
  
  inmex.meta <- readRDS("inmex_meta.rds");
  datanm.vec <- names(mdata.all)[mdata.all==1];
  dat.inx <- inmex.meta$data.lbl %in% datanm.vec;
  dat <- inmex.meta$data[, dat.inx, drop=F]; 
  
  # need to deal with missing values 
  dat <- na.omit(dat);
  
  pca3d <- list();
  if(clustOpt == "pca"){
    if(opt == "all"){
      pca <- prcomp(t(dat), center=T, scale=T);
    }else{
      dat = dat[which(rownames(dat) %in% loadEntrez),]
      pca <- prcomp(t(dat), center=T, scale=T);
    }
    
    imp.pca<-summary(pca)$importance;
    pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
    coords <- data.frame(t(signif(pca$x[,1:3], 5)));
  }else{
    require('Rtsne');
    ndat <- as.matrix(t(dat));
    max.perx <- floor((nrow(ndat)-1)/3);
    if(max.perx > 30){
      max.perx <- 30;
    }
    res <- Rtsne(ndat, dims = 3, perplexity=max.perx);
    pca3d$score$axis <- paste("t-SNE dim ", 1:3, sep="");
    coords <- data.frame(t(signif(res$Y, 5)));
  }
  
  colnames(coords) <- NULL; 
  pca3d$score$xyz <- coords;
  pca3d$score$name <- colnames(dat);
  
  facA <- as.character(inmex.meta$cls.lbl[dat.inx]);
  if(all.numeric(facA)){
    facA <- paste("Group", facA);
  }
  pca3d$score$facA <- facA;
  
  facB <-  as.character(inmex.meta$data.lbl[dat.inx]);
  if(all.numeric(facB)){
    facB <- paste("Group", facB);
  }
  pca3d$score$facB <- facB;
  
  # now set color for each group
  cols <- unique(GetColorSchema(facB));
  rgbcols <- col2rgb(cols);
  cols <- apply(rgbcols, 2, function(x){paste("rgba(", paste(x, collapse=","), ",1)", sep="")});
  pca3d$score$colors <- cols;
  
  # add shape sphere, triangles, square, pentagon (first two)
  pca3d$score$shapes <- c("sphere", "triangle");
  
  mypos <- t(coords);
  colnames(mypos) <- paste("Dim", 1:3, sep="");
  coords <- data.frame(Class=facA, Data=facB, mypos);
  write.csv(coords, file="networkanalyst_3d_pos.csv");
  
  require(RJSONIO);
  json.mat <- toJSON(pca3d, .na='null');
  sink(fileName);
  cat(json.mat);
  sink();
  current.msg <<- "Annotated data is now ready for 3D visualization!";
  return(1);
}


SaveMetaClusterLoadingJSON <- function(fileName, clustOpt, nb){
  
  inmex.meta <- readRDS("inmex_meta.rds");
  datanm.vec <- names(mdata.all)[mdata.all==1];
  nb = as.numeric(nb)
  dat.inx <- inmex.meta$data.lbl %in% datanm.vec;
  dat <- inmex.meta$data[, dat.inx, drop=F]; 
  
  # need to deal with missing values 
  dat <- na.omit(dat);
  variances = apply(dat,1, function(x){var(x)})
  df = data.frame(var = variances, inx = seq.int(1,length(variances)))
  df = df[order(-df$var),]
  inx = df$inx[c(1:nb)]
  dat = dat[inx,];
  
  pca3d <- list();
  
  pca <- prcomp(t(dat), center=T, scale=T);    
  imp.pca<-summary(pca)$importance;
  pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
  coords <- data.frame(t(signif(pca$rotation[,1:3], 5)));
  
  colnames(coords) <- NULL; 
  pca3d$score$xyz <- coords;
  pca3d$score$name <- doEntrez2SymbolMapping(rownames(pca$rotation));
  pca3d$score$entrez <- rownames(pca$rotation);
  
  loadEntrez <<- pca3d$score$entrez
  mypos <- t(coords);
  colnames(mypos) <- paste("Dim", 1:3, sep="");
  
  coords <- data.frame(mypos);
  write.csv(coords, file="networkanalyst_loadings_3d_pos.csv");
  
  require(RJSONIO);
  json.mat <- toJSON(pca3d, .na='null');
  sink(fileName);
  cat(json.mat);
  sink();
  current.msg <<- "Annotated data is now ready for 3D visualization!";
  return(1);
}

GetDatasetNamesString <- function(){
  inmex.meta <- readRDS("inmex_meta.rds");
  paste(unique(inmex.meta$data.lbl), collapse="||");
}

# this is a transient R environment for specific function, then close
PerformBatchCorrection <- function(method){
  
  print("performing combat batch correction ....");
  
  library(RSclient);
  rsc <- RS.connect();
  
  RS.assign(rsc, "my.dir", getwd()); 
  RS.eval(rsc, setwd(my.dir));
  
  my.fun <- function(){
    library('sva');
    inmex.meta <- readRDS("inmex_meta.rds");
    data.lbl <- inmex.meta$data.lbl
    pheno <- data.frame(cbind(inmex.meta$cls.lbl, data.lbl));
    modcombat <- model.matrix(~1, data=pheno)
    batch <- data.lbl;
    inmex.meta$data <- ComBat(dat=inmex.meta$data, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE);
    saveRDS(inmex.meta, "inmex_meta.rds");
    return(1);
  }
  RS.assign(rsc, my.fun);
  res <-  RS.eval(rsc, my.fun());
  RS.close(rsc);
  
  return(res);
}

##################################################
## R script for NetworkAnalyst
## Description: Functions for heatmaps
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

PerformVolcanoEnrichment<-function(file.nm, fun.type, IDs, type, inx){
  inx = as.numeric(inx)
  if(anal.type == "onedata"){
    if(dataSet$type == "array"){
      sigmat = dataSet$sig.mat
    } else {
      sigmat = dataSet$sig.mat
    }
  }else{
    sigmat = inmex.ind[selDataNm][[1]][which(inmex.ind[selDataNm][[1]][,'Pval'] < as.numeric(pvalu)),];
  }
  sigm <<- sigmat
  
  if(type == "focus"){
    gene.vec <- unlist(strsplit(IDs, "; "));
  }else if(type == "all"){
    gene.vecup <- rownames(sigmat[which(sigmat[,inx] > fcthreshu),]);
    gene.vecdown <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu),]);
    gene.vec <- c(gene.vecup, gene.vecdown);
  }else if(type == "up"){
    gene.vec <- rownames(sigmat[which(sigmat[,inx] > fcthreshu),]);
  }else{
    gene.vec <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu),]);
  }
  sym.vec <- doEntrez2SymbolMapping(gene.vec);
  names(gene.vec) <- sym.vec;
  res <- PerformEnrichAnalysis(file.nm, fun.type, gene.vec);
  return(res);
}


PerformGSEA<- function(file.nm, fun.type, netNm, mType, selectedFactorInx){
  
  
  if(tolower(fun.type) == 'kegg'){ 
    current.geneset <- LoadKEGGLib();
  }else if(tolower(fun.type) == 'reactome'){ 
    current.geneset <- LoadREACTOMELib();
  }else if(tolower(fun.type) == 'motif'){ 
    current.geneset <- LoadMotifLib();
  }else{ # GO
    current.geneset <- LoadGOLib(fun.type);
  }
  
  require("fgsea");
  
  if(anal.type == "onedata"){
    datnorm = dataSet$data.norm
    sampleNms = colnames(dataSet$data.norm);
    rankedVec<- ComputeRankedVec(dataSet, selectedFactorInx);
    
  }else{
    if(dataSet$name != selDataNm){
      dataSet <- readRDS(selDataNm);
    }
    datnorm=dataSet$data
    sampleNms = colnames(dataSet$data);
    ds = inmex.ind[selDataNm][[1]]
    rankedVec <- ComputeRankedVec(dataSet, 1);
  }
  
  fgseaRes <- fgsea(pathways = current.geneset, 
                    stats = rankedVec,
                    minSize=15,
                    maxSize=500,
                    nperm=10000);
  
  fgseaRes <- fgseaRes[!duplicated(fgseaRes$pathway),]
  
  rownames(fgseaRes) = make.names(fgseaRes$pathway, unique=TRUE)
  fgseaRes = fgseaRes[,c("size","ES", "pval", "pathway", "padj")]
  
  if(nrow(fgseaRes)<1){
    require(RJSONIO);
    SetListNms()
    initsbls <- doEntrez2SymbolMapping(list.genes);
    names(initsbls) <- list.genes
    netData <- list(sizes=listSizes, genelist=initsbls);
    netName <- paste0(netNm, ".json");
    sink(netName);
    cat(toJSON(netData));
    sink();
    return(0);
  }
  
  fgseaRes <- fgseaRes[order(fgseaRes$pval),]
  fgseaRes <<- fgseaRes
  if(nrow(fgseaRes[which(fgseaRes$pval < 0.05),])<20 ){
    if(nrow(fgseaRes)>20){
      fgseaRes = fgseaRes[c(1:20),]
    }
  }else{
    fgseaRes[which(fgseaRes$pval < 0.05),]
  } 
  
  inx <- which(names(current.geneset) %in% fgseaRes$pathway);
  current.mset = current.geneset[inx]
  current.mset = current.mset[!duplicated(names(current.mset))]
  
  ora.vec <- names(rankedVec)
  ora.nms <- doEntrez2SymbolMapping(ora.vec)
  
  hits.query <- lapply(current.mset, 
                       function(x) {
                         unique(ora.nms[ora.vec%in%unlist(x)]);
                       });
  saveRDS(hits.query, "hits_query.rds");
  set.num = unlist(lapply(current.mset, function(x){length(unique(x))}), use.names=TRUE);
  names(hits.query) <- names(current.mset);
  hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=TRUE);
  fgseaRes$hits = hit.num[which(fgseaRes$pathway  %in% names(hit.num))] 
  fgseaRes$total = set.num[which(fgseaRes$pathway %in% names(set.num))]
  
  fgseaRes = fgseaRes[which(fgseaRes$hits>1),]
  fgseaRes = fgseaRes[which(fgseaRes$hits<500),]
  fgseaRes = fgseaRes[which(fgseaRes$total<2000),]
  if(nrow(fgseaRes)<1){
    require(RJSONIO);
    SetListNms();
    initsbls = doEntrez2SymbolMapping(list.genes)
    names(initsbls) = list.genes
    netData <- list(sizes=listSizes, genelist=initsbls);
    netName = paste0(netNm, ".json");
    sink(netName);
    cat(toJSON(netData));
    sink();
    return(0);
  }
  
  fgseaRes=fgseaRes[order(fgseaRes$pval),]
  if(nrow(fgseaRes[which(fgseaRes$pval < 0.05),])<20 ){
    if(nrow(fgseaRes)>20){
      fgseaRes = fgseaRes[c(1:20),]
    }
  }else{
    fgseaRes = fgseaRes[which(fgseaRes$padj < 0.05),]
  } 
  
  fgseaRes <- data.frame(fgseaRes, stringsAsFactors=FALSE)
  
  #get gene symbols
  current.msg <<- "Functional enrichment analysis was completed";
  
  # write json
  require(RJSONIO);
  fun.anot = hits.query; 
  fun.pval = fgseaRes[,3]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  #fun.pval<-signif(fun.pval,5);  
  fun.padj = fgseaRes[,5]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  #fun.padj<-signif(fun.padj,5);  
  es.num = fgseaRes[,2]; if(length(es.num) ==1) { es.num <- matrix(es.num) };
  fun.ids <- as.vector(current.setids[names(fun.anot)]); 
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  
  json.res <- list(
    fun.link = current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    pathname = fgseaRes[,"pathway"],
    es.num = es.num,
    hits = fgseaRes[,"hits"],
    total = fgseaRes[,"total"],
    cls = dataSet$meta.info,
    sample.nms = sampleNms       
  );
  
  if(mType == "network"){
    res.mat<-matrix(0, nrow=length(fun.pval), ncol=4);
    colnames(res.mat)<-c("Total","Hits", "P.Value", "FDR");
    res.mat[,"Total"] = fgseaRes[,"total"];
    res.mat[,"Hits"] = fgseaRes[,"hits"];
    res.mat[,"P.Value"] = fgseaRes[,"pval"];
    res.mat[,"FDR"] = fgseaRes[,"padj"];
    res.mat = data.matrix(data.frame(res.mat, stringsAsFactors=FALSE));
    rownames(res.mat) = fgseaRes[,"pathway"];
    enr.mat <<- res.mat;
    list.genes <<- doEntrez2SymbolMapping(rownames(dataSet$sig.mat));
    SetListNms();
    PrepareEnrichNet(netNm, "meta", "mixed");
    file.nm <- gsub("gsea", "enrichment", file.nm)
    json.mat <- toJSON(json.res, .na='null');
    json.nm <- paste(file.nm, ".json", sep="");
  }else{
    
    json.mat <- toJSON(json.res, .na='null');
    json.nm <- paste(file.nm, ".json", sep="");
  }
  
  sink(json.nm)
  cat(json.mat);
  sink();
  
  # write csv
  
  fgseaRes <<- fgseaRes 
  ftype = fun.type
  if(fun.type %in% c("bp", "mf", "cc")){
    ftype = paste0("go_", fun.type);
  }
  csvDf <- data.frame(Name=fgseaRes$pathway, Total=fgseaRes$total, Hits=fgseaRes$hits, EnrichmentScore=fgseaRes$ES, Pval=fgseaRes$pval, Padj=fgseaRes$padj);
  write.csv(csvDf, file=paste0(file.nm, ".csv"));
  
  return(1);
}

GetQEA.pathNames<-function(){
  current.geneset <- readRDS("current_geneset.rds")
  hit.inx <- match(rownames(analSet$qea.mat),names(current.geneset));
  return(names(current.geneset)[hit.inx]);
}

PlotGSView <-function(cmpdNm,  format="png", dpi=72, width=NA){
  library("ggplot2");
  current.geneset <- readRDS("current_geneset.rds")
  nm <<- cmpdNm
  imgName <- gsub("\\/", "_",  cmpdNm);
  imgName <- gsub(" ", "_",  imgName);
  imgName <- paste(imgName, "_dpi", dpi, ".", format, sep="");
  #indx<-which(rownames(boxplot_id)==cmpdNm);
  #gene.id <- boxplot_id[indx,1];
  Cairo(file = imgName, dpi=72, width=340, height=300, type="png", bg="transparent");
  g <- plotEnrichment(current.geneset[[cmpdNm]], rankedVec)
  print(g)
  dev.off();
  return(imgName);
}

PlotGShm <-function(cmpdNm, IDs){
  ids = unlist(strsplit(IDs, "; "));
  cmpdNm <- gsub(" ", "_",  cmpdNm);
  cmpdNm <- gsub("/", "_",  cmpdNm);
  if(anal.type == "onedata"){
    subset = dataSet$data.norm[which(doEntrez2SymbolMapping(rownames(dataSet$data.norm)) %in% ids),]
    if(length(subset)<1){
      subset = dataSet$data.norm[which(rownames(dataSet$data.norm) %in% ids),]
    }
  }else{
    if(dataSet$name != selDataNm){
      dataSet <- readRDS(selDataNm);
    }
    subset = dataSet$data[which(doEntrez2SymbolMapping(rownames(dataSet$data)) %in% ids),]
    if(length(subset)<1){
      subset = dataSet$data[which(rownames(dataSet$data) %in% ids),]
    }
  }
  
  json.res <- list(
    data=subset,
    ids=doEntrez2SymbolMapping(rownames(subset)),
    entrez = rownames(subset)
  )
  library(RJSONIO);
  json.mat <- RJSONIO::toJSON(json.res, .na='null', pretty=TRUE);
  json.nm <- paste(cmpdNm,"_hm", ".json", sep="");
  
  sink(json.nm);
  cat(json.mat);
  sink();
  
  return(json.nm)
}

#For GSEA of AnalOverview page
ComputeRankedVec <- function(data, inx = 1){
  opt = rankOptGlobal;
  cls = data$cls
  if(anal.type == "metadata"){
    if(opt %in% c("mwt", "s2n", "wcx", "stu")){
      matr <- as.matrix(data$data)
    }else{
      matr <- as.matrix(readRDS("meta.resTable.rds"));
    }
  }else{
    if(opt %in% c("mwt", "s2n", "wcx", "stu")){
      matr = as.matrix(data$data.norm)
    }else{
      matr = as.matrix(data$resTable);
    }
  }
  if(opt == "mwt"){
    res = CalculateMWT(matr, cls)
    rankedVec = res$MWT
    names(rankedVec) = rownames(matr)
  }else if(opt == "stu"){
    inx1 <- which(data$cls==levels(data$cls)[1]);
    inx2 <- which(data$cls==levels(data$cls)[2]);
    res <- apply(matr, 1, function(x) {
      tmp <- try(t.test(x[inx1], x[inx2], paired = FALSE, var.equal = TRUE));
      if(class(tmp) == "try-error") {
        return(c(NA, NA));
      }else{
        return(c(tmp$statistic, tmp$p.value));
      }
    })
    res = t(res)
    rankedVec = res[,1]
    posInx = sign(rankedVec) == 1
    rankedVec[posInx] = 1000-rankedVec[posInx]
    names(rankedVec) = rownames(matr)
    
  }else if(opt == "wcx"){
    inx1 <- which(data$cls==levels(data$cls)[1]);
    inx2 <- which(data$cls==levels(data$cls)[2]);
    res <- apply(matr, 1, function(x) {
      tmp <- try(wilcox.test(x[inx1], x[inx2], paired = FALSE));
      if(class(tmp) == "try-error") {
        return(c(NA, NA));
      }else{
        return(c(tmp$statistic, tmp$p.value));
      }
    })
    res = t(res)
    rankedVec = res[,2]
    names(rankedVec) = rownames(matr)
    
  }else if (opt == "s2n"){
    res = CalculateS2N(matr, as.numeric(cls)-1)
    rankedVec = res;
  }else if(opt == "pval"){
    m = length(colnames(matr))
    if (dataSet$de.method=="limma"){
      rankedVec = as.vector(matr[,"t"]);
    } else if (dataSet$de.method=="deseq2"){
      rankedVec = as.vector(matr[,"stat"]);
    } else {
      rankedVec = as.vector(matr[,"LR"]);
    }
    names(rankedVec) = rownames(matr);
  }else{
    rankedVec = as.vector(matr[,inx]);
    names(rankedVec) = rownames(matr);
  }
  rankedVec = sort(rankedVec)
  rankedVec = rankedVec[unique(names(rankedVec))]
  rankedVec <<- rankedVec
  return(rankedVec)
}


CalculateS2N <- function(data, vec = y.vec, ...) {
  
  A <- data + 0.00000001
  
  ind1 <- which(vec==1) # cases
  n1 <- length(ind1)    
  M1 <- rowMeans(A[,ind1])
  A2 <- A*A    
  S1 <- rowMeans(A2[,ind1])   
  S1 <- S1 - M1*M1    
  S1 <- sqrt(abs((n1/(n1-1)) * S1))   
  
  ind2 <- which(vec==0) # controls
  n2 <- length(ind2)
  M2 <- rowMeans(A[,ind2])
  S2 <- rowMeans(A2[,ind2])   
  S2 <- S2 - M2*M2    
  S2 <- sqrt(abs((n2/(n2-1)) * S2))   
  
  # small sigma "fix" as used in GeneCluster
  S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
  S2 <- ifelse(S2 == 0, 0.2, S2) 
  S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
  S1 <- ifelse(S1 == 0, 0.2, S1) 
  M1 <- M1 - M2
  S1 <- S1 + S2
  s2n <- M1/S1
  
  return(s2n)
}

CalculateMWT <- function(xdat,grp,na.rm=TRUE){
  ## basic statistics
  glab = unique(grp)
  n1 = sum(grp==glab[1])
  n2 = sum(grp==glab[2])
  d1 = n1-1
  d2 = n2-1
  m1 = rowMeans(xdat[,grp==glab[1]], na.rm=na.rm)
  m2 = rowMeans(xdat[,grp==glab[2]], na.rm=na.rm)
  
  s2.g1 = rowSums((xdat[,grp==glab[1]]-m1)^2, na.rm=na.rm)/d1
  
  ## We might either have all NA in one group or variance = 0
  ## (e.g. might happen with RMA with small samples)
  ## In this situation we want to remove the gene
  s2.g1[s2.g1 == 0] <- NA
  
  s2.g2 = rowSums((xdat[,grp==glab[2]]-m2)^2, na.rm=na.rm)/d2
  s2.g2[s2.g2 == 0] <- NA
  
  ## If either s2.g1 or s2.g2 are NA this will be NA
  sig2 = (d1*s2.g1 + d2*s2.g2)/(d1+d2)
  fac = 1/n1 + 1/n2
  se2 = (sig2 * fac)
  
  ## F test
  
  lev.test = levene(xdat,grp)
  fFDR = lev.test$FDR
  fStat = lev.test$statistic
  
  
  ## ordinary Welch statistics
  se2.sep = s2.g1/n1 + s2.g2/n2
  df = se2.sep^2/((s2.g1/n1)^2/d1 + (s2.g2/n2)^2/d2)
  
  ## weighted formulas
  df.w  = fFDR*(d1+d2) + (1-fFDR)*df
  se2.w = fFDR*se2 + (1-fFDR)*se2.sep
  ds = est.hyper(z=log(se2.w),D=mean(df.w,na.rm=na.rm),d12=d1+d2)   
  
  ## ....................................... moderated Welch
  se2.com = (ds$d0*ds$s2 + df.w*se2.w)/(ds$d0 + df.w)
  Wm.t = (m1-m2)/sqrt(se2.com) ## Welch t
  df.com = ds$d0 + df.w      ## df
  
  
  Wm.pval = pt(-abs(Wm.t), df= df.com) * 2
  
  ## ................. Compute Global FDR
  
  fdr <- NULL
  
  return(list(MWT= Wm.t, coefficients=cbind((m1-m2)),pvalue = Wm.pval))
  
}

levene <- function (xdat, grp, na.rm = TRUE) 
{
  glab = unique(grp)
  ngr = length(glab)
  n = mn = s2 = NULL
  X0 = NULL
  for (i in 1:ngr) {
    ndx = grp == glab[i]
    mni = rowMeans(xdat[, ndx], na.rm = na.rm)
    x0 = xdat[, ndx] - mni
    X0 = cbind(X0, x0)
  }
  xdat = abs(X0)
  for (i in 1:ngr) {
    ndx = grp == glab[i]
    ni = sum(ndx)
    mni = rowMeans(xdat[, ndx], na.rm = na.rm)
    x0 = xdat[, ndx] - mni
    s2i = rowSums(x0 * x0, na.rm = na.rm)/(ni - 1)
    n = c(n, ni)
    mn = cbind(mn, mni)
    s2 = cbind(s2, s2i)
  }
  N = sum(n)
  mmn = rowSums(xdat, na.rm = na.rm)/N
  mn0 = mn - mmn
  num = rowSums(t(t(mn0 * mn0) * n))/(ngr - 1)
  den = rowSums(t(t(s2 * (n - 1))))/(N - ngr)
  F3 = num/den
  pval = pf(F3, df1 = ngr - 1, df2 = N - ngr)
  pval = ifelse(pval < 0.5, 2 * pval, 2 * (1 - pval))
  lvn.FDR = pval2FDR(pval)
  return(list(statistic = F3, pvalue = pval, FDR = lvn.FDR))
}

pval2FDR <-function (pval, lim = 0.7) 
{
  n1 = length(pval)
  ok.id <- 1:n1
  if (any(is.na(pval))) {
    ok.id <- which(!is.na(pval))
    pval <- na.omit(pval)
  }
  n = length(pval)
  Fp = rank(pval)/length(pval)
  p0 = sum(pval > lim)/((1 - lim) * n)
  p0 = min(p0, 1)
  FDRp = p0 * pmin(pval/Fp, 1)
  ord = order(pval)
  FDR.o = FDRp[ord]
  b = rev(cummin(rev(FDR.o)))
  FDR = rep(0, n)
  FDR[ord] = b
  out.FDR <- rep(NA, n1)
  out.FDR[ok.id] <- FDR
  attr(out.FDR, "p0") <- p0
  return(out.FDR)
}

est.hyper <- function (z, D, d12) 
{
  f2 <- function(d0, D) {
    var(z, na.rm = TRUE) - trigamma(D/2) - trigamma(d0/2)
  }
  lim = f2(100, D)
  if (lim < 0) 
    d0.est <- 100
  if (lim > 0) 
    d0.est <- uniroot(f2, c(1, 100), D = D, extendInt = "yes")$root
  s2.est <- exp(mean(z, na.rm = TRUE) - digamma(D/2) + digamma(d0.est/2) - 
                  log(d0.est/D))
  return(list(d0 = d0.est, s2 = s2.est))
}

# note: hit.query, resTable must synchronize
# ora.vec should contains entrez ids, named by their gene symbols
PerformVolcanoOneEnrichment <- function(file.nm, fun.type, IDs, inx){
  # prepare lib
  inx = as.numeric(inx);
  if(anal.type == "onedata"){
    if(dataSet$type == "array"){
      sigmat = dataSet$sig.mat
    } else {
      sigmat = dataSet$sig.mat
    }
  }else{
    sigmat = inmex.ind[selDataNm][[1]][which(inmex.ind[selDataNm][[1]][,'Pval'] < as.numeric(pvalu)),];
  }
  sigm <<- sigmat
  
  one.path.vec <- unlist(strsplit(IDs, "; "));
  
  gene.vecup <- rownames(sigmat[which(sigmat[,inx] > fcthreshu),]);
  gene.vecdown <- rownames(sigmat[which(sigmat[,inx] < -fcthreshu),]);
  ora.vec <- c(gene.vecup, gene.vecdown);
  
  
  sym.vec <- doEntrez2SymbolMapping(ora.vec);
  names(ora.vec) <- sym.vec;
  
  current.geneset <- list()
  current.geneset[["Set"]] = one.path.vec
  current.geneset[["Set2"]] = one.path.vec
  
  # prepare query
  ora.nms <- names(ora.vec);
  
  # prepare for the result table
  set.size<-length(current.geneset);
  res.mat<-matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat)<-names(current.geneset);
  colnames(res.mat)<-c("Total", "Expected", "Hits", "P.Value", "FDR");
  
  # need to cut to the universe covered by the pathways, not all genes
  current.universe <- unique(unlist(current.geneset)); 
  hits.inx <- ora.vec %in% current.universe;
  ora.vec <- ora.vec[hits.inx];
  ora.nms <- ora.nms[hits.inx];
  
  q.size<-length(ora.vec);
  
  # get the matched query for each pathway
  hits.query <- lapply(current.geneset, 
                       function(x) {
                         ora.nms[ora.vec%in%unlist(x)];
                       }
  );
  
  saveRDS(hits.query, "hits_query.rds");
  
  names(hits.query) <- names(current.geneset);
  hit.num<-unlist(lapply(hits.query, function(x){length(unique(x))}), use.names=FALSE);
  
  # total unique gene number
  #uniq.count <- length(current.universe);
  uniq.count <- nrow(dataSet$data.norm);
  
  # unique gene count in each pathway
  set.size <- unlist(lapply(current.geneset, length));
  
  res.mat[,1]<-set.size;
  res.mat[,2]<-q.size*(set.size/uniq.count);
  res.mat[,3]<-hit.num;
  
  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
  res.mat[,4]<- raw.pvals;
  res.mat[,5] <- raw.pvals;
  
  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];
  if(nrow(res.mat)> 0){
    # order by p value
    ord.inx<-order(res.mat[,4]);
    #res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];
    
    imp.inx <- res.mat[,4] <= 0.05;
    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      if(sum(imp.inx) > 120){
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }else{
    return(0);
  }
  
  #get gene symbols
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
  res.mat[,"Hits"] = res.mat[,"Hits"]
  enr.mat <<- res.mat
  current.msg <<- "Functional enrichment analysis was completed";
  
  # write json
  require(RJSONIO);
  fun.anot = hits.query; 
  total = resTable[,2]; if(length(total) ==1) { total <- matrix(total) };
  fun.pval = resTable[,5]; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  fun.padj = resTable[,6]; if(length(fun.padj) ==1) { fun.padj <- matrix(fun.padj) };
  hit.num = resTable[,4]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
  fun.ids <- as.vector(current.setids[names(fun.anot)]); 
  if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  json.res <- list(
    fun.link = current.setlink[1],
    fun.anot = fun.anot,
    fun.ids = fun.ids,
    fun.pval = fun.pval,
    fun.padj = fun.padj,
    hit.num = hit.num,
    total= total
  );
  json.mat <- toJSON(json.res, .na='null');
  json.nm <- paste(file.nm, ".json", sep="");
  
  sink(json.nm)
  cat(json.mat);
  sink();
  
  # write csv
  fun.hits <<- hits.query;
  fun.pval <<- resTable[,5];
  hit.num <<- resTable[,4];
  csv.nm <- paste(file.nm, ".csv", sep="");    
  write.csv(resTable, file=csv.nm, row.names=F);
  
  return(1);
}

##################################################
## R script for NetworkAnalyst
## Description: General graph manipulation functions 
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################

DecomposeGraph <- function(gObj, minNodeNum = 3){
  # now decompose to individual connected subnetworks
  comps <-decompose.graph(gObj, min.vertices=minNodeNum);
  if(length(comps) == 0){
    current.msg <<- paste("No subnetwork was identified with at least", minNodeNum, "nodes!");
    return(NULL);
  }
  
  # first compute subnet stats
  net.stats <- ComputeSubnetStats(comps);
  ord.inx <- order(net.stats[,1], decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  comps <- comps[ord.inx];
  names(comps) <- rownames(net.stats) <- paste("subnetwork", 1:length(comps), sep="");
  
  # note, we report stats for all nets (at least 3 nodes);
  hit.inx <- net.stats$Node >= minNodeNum;
  comps <- comps[hit.inx];
  
  # now record
  ppi.comps <<- comps;
  net.stats <<- net.stats;
  
  sub.stats <- unlist(lapply(comps, vcount)); 
  return(sub.stats);
}


ComputeSubnetStats <- function(comps){
  net.stats <- as.data.frame(matrix(0, ncol = 3, nrow = length(comps)));
  colnames(net.stats) <- c("Node", "Edge", "Query");
  for(i in 1:length(comps)){
    g <- comps[[i]];
    net.stats[i,] <- c(vcount(g),ecount(g),sum(seed.proteins %in% V(g)$name));
  }
  return(net.stats);
}


# new range [a, b]
rescale2NewRange <- function(qvec, a, b){
  q.min <- min(qvec);
  q.max <- max(qvec);
  if(length(qvec) < 50){
    a <- a*2;
  }
  if(q.max == q.min){
    new.vec <- rep(8, length(qvec));
  }else{
    coef.a <- (b-a)/(q.max-q.min);
    const.b <- b - coef.a*q.max;
    new.vec <- coef.a*qvec + const.b;
  }
  return(new.vec);
}


ComputeColorGradient <- function(nd.vec, background="black", centered, colorblind){
  require("RColorBrewer");
  
  minval = min(nd.vec, na.rm=TRUE);
  maxval = max(nd.vec, na.rm=TRUE);
  res = maxval-minval;
  
  if(res == 0){
    return(rep("#FF0000", length(nd.vec)));
  }
  color <- GetColorGradient(background, centered, colorblind);
  breaks <- generate_breaks(nd.vec, length(color), center = centered);
  return(scale_vec_colours(nd.vec, col = color, breaks = breaks));
}

generate_breaks = function(x, n, center = F){
  if(center){
    m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    res = seq(-m, m, length.out = n + 1)
  }
  else{
    res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
  }
  return(res)
}

scale_vec_colours = function(x, col = rainbow(10), breaks = NA){
  breaks <- sort(unique(breaks));
  return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
}

GetColorGradient <- function(background, center, colorblind=F) {
  if (background == "black") {
    if (center) {
      if (colorblind) {
        return(c(colorRampPalette(c("#3182bd", "#6baed6", "#bdd7e7", "#eff3ff"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      } else {
        return(c(colorRampPalette(c("#31A231", "#5BC85B", "#90EE90", "#C1FFC1"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      }
    } else {
      if (colorblind) {
        return(c(colorRampPalette(c("#3182bd", "#bbfdff"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      } else {
        return(colorRampPalette(rev(heat.colors(9)))(100))
      }
    }
  } else {
    if (center) {
      if (colorblind) {
        return(c(colorRampPalette(c("#3182bd", "#bbfdff"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)))
      } else {
        return(c(colorRampPalette(c("#137B13", "#31A231", "#5BC85B", "#90EE90"))(50), colorRampPalette(c("#FF7783", "#E32636", "#BD0313", "#96000D"))(50)))
      }
    } else {
      return(colorRampPalette(hsv(h = seq(0.72, 1, 0.035), s = 0.72, v = 1))(100))
    }
  }
}

# also save to GraphML
ExportNetwork <- function(fileName){
  current.net <- ppi.comps[[current.net.nm]];
  write.graph(current.net, file=fileName, format="graphml");
}


ExtractModule<- function(nodeids){
  set.seed(8574);
  nodes <- strsplit(nodeids, ";")[[1]];
  
  g <- ppi.comps[[current.net.nm]];
  # try to see if the nodes themselves are already connected
  hit.inx <- V(g)$name %in% nodes; 
  gObj <- induced.subgraph(g, V(g)$name[hit.inx]);
  
  # now find connected components
  comps <-decompose.graph(gObj, min.vertices=1);
  
  if(length(comps) == 1){ # nodes are all connected
    g <- comps[[1]];
  }else{
    # extract modules
    paths.list <-list();
    sd.len <- length(nodes);
    for(pos in 1:sd.len){
      paths.list[[pos]] <- get.shortest.paths(g, nodes[pos], nodes[-(1:pos)])$vpath;
    }
    nds.inxs <- unique(unlist(paths.list));
    nodes2rm <- V(g)$name[-nds.inxs];
    g <- simplify(delete.vertices(g, nodes2rm));
  }
  nodeList <- get.data.frame(g, "vertices");
  if(nrow(nodeList) < 3){
    return ("NA");
  }
  
  module.count <- module.count + 1;
  module.nm <- paste("module", module.count, sep="");
  colnames(nodeList) <- c("Id", "Label");
  ndFileNm = paste(module.nm, "_node_list.csv", sep="");
  write.csv(nodeList, file=ndFileNm, row.names=F, quote=F);
  
  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  edgFileNm = paste(module.nm, "_edge_list.csv", sep="");
  write.csv(edgeList, file=edgFileNm, row.names=F, quote=F);
  
  filenm <- paste(module.nm, ".json", sep="");
  
  # record the module 
  ppi.comps[[module.nm]] <<- g;
  UpdateSubnetStats();
  
  module.count <<- module.count
  
  convertIgraph2JSON(module.nm, filenm);
  return (filenm);
}

PerformLayOut <- function(net.nm, algo){
  g <- ppi.comps[[net.nm]];
  vc <- vcount(g);
  if(algo == "Default"){
    if(vc > 3000) {
      pos.xy <- layout.lgl(g, maxiter = 100);
    }else if(vc > 2000) {
      pos.xy <- layout.lgl(g, maxiter = 150);
    }else if(vc > 1000) {
      pos.xy <- layout.lgl(g, maxiter = 200);
    }else if(vc < 150){
      pos.xy <- layout.kamada.kawai(g);
    }else{
      pos.xy <- layout.fruchterman.reingold(g);
    }
  }else if(algo == "FrR"){
    pos.xy <- layout.fruchterman.reingold(g);
  }else if(algo == "random"){
    pos.xy <- layout.random(g);
  }else if(algo == "lgl"){
    if(vc > 3000) {
      pos.xy <- layout.lgl(g, maxiter = 100);
    }else if(vc > 2000) {
      pos.xy <- layout.lgl(g, maxiter = 150);
    }else {
      pos.xy <- layout.lgl(g, maxiter = 200);
    }
  }else if(algo == "gopt"){
    # this is a slow one
    if(vc > 3000) {
      maxiter = 50;
    }else if(vc > 2000) {
      maxiter = 100;
    }else if(vc > 1000) {
      maxiter = 200;
    }else{
      maxiter = 500;
    }
    pos.xy <- layout.graphopt(g, niter=maxiter);
  }
  pos.xy;
}

UpdateNetworkLayout <- function(algo, filenm){
  current.net <- ppi.comps[[current.net.nm]];
  #pos.xy <- PerformLayOut_mem(current.net.nm, algo);
  pos.xy <- PerformLayOut(current.net.nm, algo);
  nms <- V(current.net)$name;
  nodes <- vector(mode="list");
  for(i in 1:length(nms)){
    nodes[[i]] <- list(
      id=nms[i], 
      x=pos.xy[i,1], 
      y=pos.xy[i,2]
    );
  }
  # now only save the node pos to json
  require(RJSONIO);
  netData <- list(nodes=nodes);
  sink(filenm);
  cat(toJSON(netData));
  sink();
  return(filenm);
}

getGraphStatsFromFile <- function(){
  g <- ppi.comps[[net.nmu]];
  nms <- V(g)$name;
  edge.mat <- get.edgelist(g);
  return(c(length(nms), nrow(edge.mat)));        
}

GetNetUploaded <- function(){
  netUploadU
}

GetNetsNameString <- function(){
  paste(rownames(net.stats), collapse="||");
}

# support walktrap, infomap and lab propagation
FindCommunities <- function(method="walktrap", use.weight=FALSE){
  
  # make sure this is the connected
  current.net <- ppi.comps[[current.net.nm]];
  g <- current.net;
  if(!is.connected(g)){
    g <- decompose.graph(current.net, min.vertices=2)[[1]];
  }
  total.size <- length(V(g));
  
  if(use.weight){ # this is only tested for walktrap, should work for other method
    # now need to compute weights for edges
    egs <- get.edges(g, E(g)); #node inx
    nodes <- V(g)$name;
    # conver to node id
    negs <- cbind(nodes[egs[,1]],nodes[egs[,2]]);
    
    # get min FC change
    base.wt <- min(abs(seed.expr))/10;
    
    # check if user only give a gene list without logFC or all same fake value
    if(length(unique(seed.expr)) == 1){
      seed.expr <- rep(1, nrow(negs));
      base.wt <- 0.1; # weight cannot be 0 in walktrap
    }
    
    wts <- matrix(base.wt, ncol=2, nrow = nrow(negs));
    for(i in 1:ncol(negs)){
      nd.ids <- negs[,i];
      hit.inx <- match(names(seed.expr), nd.ids);
      pos.inx <- hit.inx[!is.na(hit.inx)];
      wts[pos.inx,i]<- seed.expr[!is.na(hit.inx)]+0.1;
    }
    nwt <- apply(abs(wts), 1, function(x){mean(x)^2})    
  }
  
  if(method == "walktrap"){
    fc <- walktrap.community(g);
  }else if(method == "infomap"){
    fc <- infomap.community(g);
  }else if(method == "labelprop"){
    fc <- label.propagation.community(g);
  }else{
    print(paste("Unknown method:", method));
    return ("NA||Unknown method!");
  }
  
  if(length(fc) == 0 || modularity(fc) == 0){
    return ("NA||No communities were detected!");
  }
  
  # only get communities
  communities <- communities(fc);
  community.vec <- vector(mode="character", length=length(communities));
  gene.community <- NULL;
  qnum.vec <- NULL;
  pval.vec <- NULL;
  rowcount <- 0;
  nms <- V(g)$name;
  hit.inx <- match(nms, ppi.net$node.data[,1]);
  sybls <- ppi.net$node.data[hit.inx,2];
  names(sybls) <- V(g)$name;
  for(i in 1:length(communities)){
    # update for igraph 1.0.1 
    path.ids <- communities[[i]];
    psize <- length(path.ids);
    if(psize < 5){
      next; # ignore very small community
    }
    hits <- seed.proteins %in% path.ids;
    qnums <- sum(hits);
    if(qnums == 0){
      next; # ignor community containing no queries
    }
    
    rowcount <- rowcount + 1;
    pids <- paste(path.ids, collapse="->");
    #path.sybls <- V(g)$Label[path.inx];
    path.sybls <- sybls[path.ids];
    com.mat <- cbind(path.ids, path.sybls, rep(i, length(path.ids)));
    gene.community <- rbind(gene.community, com.mat);
    qnum.vec <- c(qnum.vec, qnums);
    
    # calculate p values (comparing in- out- degrees)
    #subgraph <- induced.subgraph(g, path.inx);
    subgraph <- induced.subgraph(g, path.ids);
    in.degrees <- degree(subgraph);
    #out.degrees <- degree(g, path.inx) - in.degrees;
    out.degrees <- degree(g, path.ids) - in.degrees;
    ppval <- wilcox.test(in.degrees, out.degrees)$p.value;
    ppval <- signif(ppval, 3);
    pval.vec <- c(pval.vec, ppval);
    
    # calculate community score
    community.vec[rowcount] <- paste(c(psize, qnums, ppval, pids), collapse=";");
  }
  
  ord.inx <- order(pval.vec, decreasing=F);
  community.vec <- community.vec[ord.inx];
  qnum.vec <- qnum.vec[ord.inx];
  ord.inx <- order(qnum.vec, decreasing=T);
  community.vec <- community.vec[ord.inx];
  
  all.communities <- paste(community.vec, collapse="||");
  colnames(gene.community) <- c("Id", "Label", "Module");
  write.csv(gene.community, file="module_table.csv", row.names=F);
  return(all.communities);
}

community.significance.test <- function(graph, vs, ...) {
  subgraph <- induced.subgraph(graph, vs)
  in.degrees <- degree(subgraph)
  out.degrees <- degree(graph, vs) - in.degrees
  wilcox.test(in.degrees, out.degrees, ...)
}


convertIgraph2JSON <- function(net.nm, filenm){
  g <- ppi.comps[[net.nm]];
  # annotation
  nms <- V(g)$name;
  hit.inx <- match(nms, ppi.net$node.data[,1]);
  lbls <- ppi.net$node.data[hit.inx,2];
  
  # setup shape (gene circle, other squares)
  shapes <- rep("circle", length(nms));
  
  # get edge data
  edge.mat <- get.edgelist(g);
  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2]);
  
  # now get coords
  #pos.xy <- PerformLayOut_mem(net.nm, "Default");
  pos.xy <- PerformLayOut(net.nm, "Default");
  # get the note data
  node.btw <- as.numeric(betweenness(g));
  node.dgr <- as.numeric(degree(g));
  node.exp <- as.vector(expr.vec[nms]);
  
  node.exp[is.na(node.exp)] = 0;
  # node size to degree values
  if(vcount(g)>500){
    min.size = 1;
  }else if(vcount(g)>200){
    min.size = 2;
  }else{
    min.size = 3;
  }
  node.sizes <- as.numeric(rescale2NewRange((log(node.dgr))^2, min.size, 9));
  centered = T;
  notcentered = F;
  # update node color based on betweenness
  require("RColorBrewer");
  topo.val <- log(node.btw+1);
  topo.colsb <- ComputeColorGradient(topo.val, "black", notcentered, FALSE);
  topo.colsw <-  ComputeColorGradient(topo.val, "white", notcentered, FALSE);
  topo.colsc <-  ComputeColorGradient(topo.val, "colorblind", notcentered, TRUE);
  
  # color based on expression
  bad.inx <- is.na(node.exp) | node.exp==0;
  if(!all(bad.inx)){
    exp.val <- node.exp;
    node.colsb.exp <- ComputeColorGradient(exp.val, "black", centered, FALSE); 
    node.colsw.exp <- ComputeColorGradient(exp.val, "white", centered, FALSE);
    node.colsc.exp <- ComputeColorGradient(exp.val, "colorblind", centered, TRUE);
    
    node.colsb.exp[bad.inx] <- "#d3d3d3"; 
    node.colsw.exp[bad.inx] <- "#c6c6c6"; 
    node.colsc.exp[bad.inx] <- "#c6c6c6"; 
    # node.colsw.exp[bad.inx] <- "#b3b3b3";
  }else{
    node.colsb.exp <- rep("#d3d3d3",length(node.exp)); 
    node.colsw.exp <- rep("#c6c6c6",length(node.exp)); 
    node.colsc.exp <- rep("#b3b3b3",length(node.exp)); 
  }
  
  # now update for bipartite network
  notbipartite = c("ppi", "tissuecoex", "cellcoex", "tissueppi")
  if(!ppi.net$db.type %in% c("ppi", "tissuecoex", "cellcoex", "tissueppi")){ # the other part miRNA or TF will be in square
    if(ppi.net$db.typ == "tfmir"){
      db.path <- paste(lib.path, data.org, "/mirlist.rds", sep="");
      db.map <-  readRDS(db.path);
      mir.inx <- nms %in% db.map[,1];
      tf.inx <- nms %in% edge.mat[,"source"];
      shapes[mir.inx] <- "square";
      shapes[tf.inx] <- "diamond";
      node.sizes[mir.inx] <- node.sizes[mir.inx] + 0.5;
      
      # update mir node color
      node.colsw.exp[tf.inx] <- topo.colsw[tf.inx] <- "#00ff00";
      node.colsw.exp[tf.inx] <- topo.colsw[tf.inx] <- "#00ff00";
      node.colsc.exp[tf.inx] <- topo.colsc[tf.inx] <- "#00ff00";
      node.colsw.exp[mir.inx] <- topo.colsw[mir.inx] <- "#306EFF"; # dark blue
      node.colsb.exp[mir.inx] <- topo.colsb[mir.inx] <- "#98F5FF";
      node.colsc.exp[mir.inx] <- topo.colsb[mir.inx] <- "#98F5FF";
    }else{
      mir.inx <- nms %in% edge.mat[,"target"];
      shapes[mir.inx] <- "square";
      node.sizes[mir.inx] <- node.sizes[mir.inx] + 0.5;
      
      # update mir node color
      node.colsw.exp[mir.inx] <- topo.colsw[mir.inx] <- "#306EFF"; # dark blue
      node.colsb.exp[mir.inx] <- topo.colsb[mir.inx] <- "#98F5FF";
      node.colsc.exp[mir.inx] <- topo.colsc[mir.inx] <- "#98F5FF";
    }
  }
  seed.inx <- nms %in% unique(seed.proteins);
  seed_arr <- rep("notSeed",length(node.dgr));
  seed_arr[seed.inx] <- "seed";
  
  # now create the json object
  nodes <- vector(mode="list");
  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id=nms[i],
      idnb = i, 
      label=lbls[i],
      x = pos.xy[i,1],
      y = pos.xy[i,2],
      size=node.sizes[i], 
      seedArr = seed_arr[i],
      type=shapes[i],
      colorb=topo.colsb[i],
      colorw=topo.colsw[i],
      colorc=topo.colsc[i],
      expcolb=node.colsb.exp[i],
      expcolw=node.colsw.exp[i],
      expcolc=node.colsc.exp[i],
      
      highlight = 0,
      attributes=list(
        expr = node.exp[i],
        degree=node.dgr[i], 
        between=node.btw[i])
    );
  }
  
  # save node table
  nd.tbl <- data.frame(Id=nms, Label=lbls, Degree=node.dgr, Betweenness=round(node.btw,2), Expression=node.exp);
  # order 
  ord.inx <- order(nd.tbl[,3], nd.tbl[,4], decreasing = TRUE)
  nd.tbl <- nd.tbl[ord.inx, ];
  write.csv(nd.tbl, file="node_table.csv", row.names=FALSE);
  
  # covert to json
  require(RJSONIO);
  netData <- list(nodes=nodes, edges=edge.mat);
  netUploadU <<-0
  sink(filenm);
  cat(toJSON(netData));
  sink();
}


# from to should be valid nodeIDs
GetShortestPaths <- function(from, to){
  current.net <- ppi.comps[[current.net.nm]];
  paths <- get.all.shortest.paths(current.net, from, to)$res;
  if(length(paths) == 0){
    return (paste("No connection between the two nodes!"));
  }
  
  path.vec <- vector(mode="character", length=length(paths));
  for(i in 1:length(paths)){
    path.inx <- paths[[i]]; 
    path.ids <- V(current.net)$name[path.inx];
    path.sybls <- path.ids;
    pids <- paste(path.ids, collapse="->");
    psbls <- paste(path.sybls, collapse="->");
    path.vec[i] <- paste(c(pids, psbls), collapse=";")
  }
  
  if(length(path.vec) > 50){
    path.vec <- path.vec[1:50];
  }
  
  all.paths <- paste(path.vec, collapse="||");
  return(all.paths);
}

GetNetsName <- function(){
  rownames(net.stats);
}

GetNetsEdgeNum <- function(){
  as.numeric(net.stats$Edge);
}

GetNetsNodeNum <- function(){
  as.numeric(net.stats$Node);
}

GetNetsQueryNum <- function(){
  as.numeric(net.stats$Query);
}

##################################################
## R script for NetworkAnalyst
## Description: Graph IO functions for network upload module
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################

ReadGraphFile <- function(fileName, fileType) {
  require("igraph");
  types_arr <<- "";
  
  fileTypeu <<- fileType;
  current.msg <<- NULL;
  
  if(fileType == "graphml"){
    graphX = tryCatch({
      read_graph(fileName, format = "graphml")
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, finally = {
      
    })
  }else if(fileType == "sif"){
    graphX = tryCatch({
      read.sif(fileName, format="igraph", directed = FALSE, sep="\t")
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, finally = {
      
    })
  }else if(fileType == "txt"){
    df <- read.table(fileName, header=FALSE, stringsAsFactors = FALSE)
    df = as.matrix(df)
    graphX = tryCatch({
      graph_from_edgelist(df)
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, finally = {
      
    })
  }else if(fileType == "json"){
    require("RJSONIO");
    dat = fromJSON(fileName);
    dfn = unlist(dat$elements$nodes);
    conv = data.frame(id1=dfn[which(names(dfn)=='data.id')], name1=dfn[which(names(dfn)=='data.name')]);
    dfe = unlist(dat$elements$edges);
    dffe = data.frame(id1=dfe[which(names(dfe) == "data.source")], id2=dfe[which(names(dfe) == "data.target")]);
    dfint = merge(conv, dffe, by="id1");
    colnames(conv) = c("id2", "name2");
    df = merge(conv, dfint, by="id2");
    df = df[,c("name1", "name2")];
    df=as.matrix(df)
    
    graphX = tryCatch({
      graph_from_edgelist(df, directed=FALSE)
    }, warning = function(w) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, error = function(e) {
      current.msg <<- "Wrong format, please make sure that the file is formatted correctly";
      return(0)
    }, finally = {
      
    })
  }else{
    current.msg <<- "Unknown format, please make sure that the file is saved in the supported formats!";
    return(0)
  }
  
  if(!is_igraph(graphX)){
    current.msg <<- "Failed to parse your file, please make sure that the file is formatted correctly";
    return(0)
  }
  current.msg <<- "Sucessfully parsed your graph file!";
  print(current.msg);
  nms <- V(graphX)$name;
  if(length(nms)<1){
    nms <- V(graphX)$id;
    graphX = set_vertex_attr(graphX, "name", value=nms)
  }
  node.data = data.frame(nms, nms);
  seed.proteins <<- nms;
  seed.genes <<- seed.proteins;
  e=get.edgelist(graphX)
  edge.data= data.frame(Source=e[,1], Target=e[,2])
  seed.expr <<- rep(0, length(node.data));
  substats <- DecomposeGraph(graphX);
  net.nm <- names(ppi.comps)[1];
  net.nmu <<- net.nm;
  current.net.nm <<- net.nm;
  g <- ppi.comps[[net.nm]];
  ppi.net <<- list(db.type="abc",
                   db.type="ppi", 
                   order=1, 
                   seeds=nms, 
                   table.nm=" ", 
                   node.data = node.data,
                   edge.data = edge.data
  );
  
  convertIgraph2JSONFromFile(net.nm, "networkanalyst_0.json");
  return(1);
}


# create igraph from the edgelist saved from graph DB
# and decompose into subnets

convertIgraph2JSONFromFile <- function(net.nm, filenm){
  
  g <- ppi.comps[[net.nm]];
  
  # annotation
  nms <- V(g)$name;
  
  # setup shape (gene circle, other squares)
  shapes <- rep("circle", length(nms));
  
  # get edge data
  edge.mat <- get.edgelist(g);
  
  edge.mat1 = data.frame(edge.mat)
  edge.mat1$color = "target"
  edge.mat1 = as.matrix(edge.mat1)
  
  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], color = edge.mat1[,3]);
  
  # get the note data
  node.btw <- as.numeric(betweenness(g));
  #node.clo <- as.numeric(closeness(g));
  #node.adh <- as.numeric(adhesion(g));
  node.eig <- eigen_centrality(g);
  node.eig = as.numeric(node.eig$vector);
  node.tra=transitivity(g,type=c("local"))
  
  node.dgr <- as.numeric(degree(g));
  node.exp <- as.numeric(get.vertex.attribute(g, name="Expression", index = V(g)));
  
  if(length(node.exp) == 0){
    node.exp <- rep(0,length(node.dgr)); 
  }
  
  # node size to degree values
  if(vcount(g)>500){
    min.size = 1;
  }else if(vcount(g)>200){
    min.size = 2;
  }else{
    min.size = 3;
  }
  
  minval = min(node.dgr, na.rm=T);
  maxval = max(node.dgr, na.rm=T);
  result = maxval-minval;
  
  if(result == 0){
    node.sizes <- rep((log(node.dgr))^2, length(nms));
  }else{
    node.sizes <- as.numeric(rescale2NewRange((log(node.dgr))^2, min.size, 9));
  }
  
  centered = T;
  notcentered = F;
  
  # update node color based on betweenness
  require("RColorBrewer");
  topo.val <- log(node.btw+1);
  topo.colsb <- ComputeColorGradient(topo.val, "black", notcentered, FALSE);
  topo.colsw <-  ComputeColorGradient(topo.val, "white", notcentered, FALSE);
  topo.colsc <-  ComputeColorGradient(topo.val, "colorblind", notcentered, TRUE);
  
  # color based on expression
  bad.inx <- is.na(node.exp) | node.exp==0;
  if(!all(bad.inx)){
    exp.val <- node.exp;
    node.colsb.exp <- ComputeColorGradient(exp.val, "black", centered, FALSE); 
    node.colsw.exp <- ComputeColorGradient(exp.val, "white", centered, FALSE);
    node.colsc.exp <- ComputeColorGradient(exp.val, "colorblind", centered, TRUE);
    node.colsb.exp[bad.inx] <- "#d3d3d3"; 
    node.colsw.exp[bad.inx] <- "#c6c6c6"; 
    node.colsc.exp[bad.inx] <- "#99ddff";
  }else{
    node.colsb.exp <- rep("#d3d3d3",length(node.exp)); 
    node.colsw.exp <- rep("#c6c6c6",length(node.exp)); 
    node.colsc.exp <- rep("#99ddff",length(node.exp)); 
  }
  
  # now update for bipartite network
  gene.inx <- nms %in% edge.mat[,"source"];
  mir.inx <- nms %in% edge.mat[,"target"];
  node_attr = list.vertex.attributes(g);
  
  attr=list();
  for(j in 1:length(node_attr)){
    attr[[node_attr[j]]] = vertex_attr(g, node_attr[j])
  }
  attr_names = names(attr);
  attr_nd = list();
  arr = list()
  for(i in 1:length(node.sizes)){
    for(j in 1:length(attr)){
      attr_nd[node_attr[j]] = as.character(unlist(attr[node_attr[j]])[i])
    }
    arr[[i]] = attr_nd;
  }
  network_prop = list();
  for(i in 1:length(node.sizes)){
    network_prop[[i]]  <- list(
      #Closeness = node.clo[i],
      Eigen = node.eig[i],
      Transitivity = node.tra[i]
    )
  }
  pos.xy <- PerformLayOut(net.nm, "Default");
  lblsu <<- nms;
  # now create the json object
  nodes <- vector(mode="list");
  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id=nms[i],
      idnb=i,
      x = pos.xy[i,1],
      y = pos.xy[i,2],
      label=nms[i],
      size=node.sizes[i], 
      type="circle",
      colorb=topo.colsb[i],
      colorw=topo.colsw[i],
      colorc=topo.colsc[i],
      topocolb=topo.colsb[i],
      topocolw=topo.colsw[i],
      topocolc=topo.colsc[i],
      expcolb=node.colsb.exp[i],
      expcolw=node.colsw.exp[i],
      expcolc=node.colsc.exp[i],
      user=c(arr[[i]], network_prop[[i]]),
      attributes=list(
        expr = node.exp[i],
        degree=node.dgr[i], 
        between=node.btw[i]
      )
    );
  }
  
  # save node table
  nd.tbl <- data.frame(Id=nms, Degree=node.dgr, Betweenness=round(node.btw,2));
  # order 
  ord.inx <- order(nd.tbl[,2], nd.tbl[,3], decreasing = TRUE)
  nd.tbl <- nd.tbl[ord.inx, ];
  write.csv(nd.tbl, file="node_table.csv", row.names=FALSE);
  # covert to json
  require(RJSONIO);
  netData <- list(nodes=nodes, edges=edge.mat);
  netUploadU <<-1
  sink(filenm);
  cat(toJSON(netData));
  sink();
}


read.sif <- function (sif.file, format = "graphNEL", directed = FALSE, header = TRUE, sep = "\t", ...) {
  
  net <- read.csv(file = sif.file, sep = sep, colClasses = "character", header = header, ...)
  
  # Assume form: node1 linktype node2 side.info..
  if ( ncol(net) > 2 ) { 
    
    # remove NA nodes 
    nas <- apply(net, 1, function (x) {any(is.na(x[c(1,3)]))})
    if (any(nas)) {
      net <- net[!nas, ]
      warning("NAs removed from network node list, ", sum(nas), " edges removed.")
    }
    
    net <- graph.edgelist(as.matrix(net[, -2]), directed = directed)
    
  } else if ( ncol(net) == 2 ) { # assume form: node1 node2
    
    # remove NA nodes 
    nas <- apply(net, 1, function (x) {any(is.na(x))})
    if (any(nas)) {
      net <- net[!nas, ]
      warning("NAs removed from network node list, ", sum(nas), " edges removed.")
    }
    
    net <- graph.edgelist(cbind(net[,1],net[,2]), directed = directed)
  }
  
  if (format == "graphNEL") { net <- igraph.to.graphNEL(net) }
  
  return(net);
}

##################################################
## R scripts for NetworkAnalyst 
## Description: Meta Analysis Methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# for multiple class, only select two
# also record all grp lbls
SetGroupContrast <- function(dataName, grps){
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  if(length(levels(dataSet$cls))>2){ 
    print("Updating group contrasts .....");
    grp.nms <- strsplit(grps, " vs. ")[[1]];
    sel.inx <- as.character(dataSet$cls) %in% grp.nms;
    
    # regenerate factor to drop levels, force the levels order
    group <- factor(dataSet$cls[sel.inx], levels=grp.nms);  
    data <- dataSet$data[, sel.inx];
    dataSet$cls <- group;
    dataSet$data <- data;
    RegisterData(dataSet);  
  }
}

# determine if all annotated data are ready for meta-analysis
CheckMetaDataIntegrity<-function(){
  
  performedDE <<- FALSE;
  
  if(length(mdata.all) == 0){
    current.msg <<-"Please upload your data or try our example datasets!";
    print(current.msg);
    return(0);
  }
  
  include.inx <- mdata.all==1;
  if(sum(include.inx) == 0){
    current.msg <<-"No dataset is selected for analysis!";
    print(current.msg);
    return(0);
  }
  
  sel.nms <- names(mdata.all)[include.inx];
  clss <- list();
  if(meta.upload){
    # update meta data only for select/deselect datasets
    inmex.meta.orig <- readRDS("inmex.meta.orig.rds");
    hit.inx <- inmex.meta.orig$data.lbl %in% sel.nms;
    data <- inmex.meta.orig$data[, hit.inx];
    id.type <- inmex.meta.orig$id.type;
    cls.lbl <- factor(inmex.meta.orig$cls.lbl[hit.inx]);
    data.lbl <- factor(inmex.meta.orig$data.lbl[hit.inx]);
    common.matrix <- data;
    nms.vec <<- rownames(inmex.meta.orig$data);
    smps.vec <<- data.lbl
  }else{   
    # first check that all class labels are consistent
    dataSet <- readRDS(sel.nms[1]);
    lvls <- levels(dataSet$cls);
    id.type <- dataSet$id.type;
    clss[[1]] = dataSet$cls;
    nms <- rownames(dataSet$data);
    shared.nms <- nms;
    for(i in 2:length(sel.nms)){
      dataSet <- readRDS(sel.nms[i]);
      clss[[i]] = dataSet$cls;
      # check if class label is consistent
      if(!all(levels(dataSet$cls) == lvls)){
        current.msg <<- paste(sel.nms[i], "has different group labels", 
                              paste(levels(dataSet$cls), collapse=":"), 
                              "from", sel.nms[1], paste(lvls, collapse=":"));
        print(current.msg);
        return(0);
      }
      
      # check and record if there is common genes            
      shared.nms <- intersect(shared.nms, rownames(dataSet$data));
      if(length(shared.nms) < 10){
        current.msg <<- paste(sel.nms[i], "has less than 10 common genes/probes from previous data sets");
        print(current.msg);
        return(0);
      }
      
      # check gene id type
      if(dataSet$id.type != id.type){
        current.msg <<- paste(sel.nms[i], "has different gene/probe ID from", sel.nms[1]);
        print(current.msg);
        return(0);
      }
    }
    
    nrepu<<-unlist(lapply(clss,FUN=function(x)length(x)));      
    
    print("Passed exp condition check!");
    
    # now construct a common matrix to faciliated plotting across all studies
    dataName <- sel.nms[1];
    dataSet <- readRDS(dataName);
    common.matrix <- dataSet$data[shared.nms, ];
    nms.vec = rownames(dataSet$data);
    smps.vec = colnames(dataSet$data);
    data.lbl <- rep(dataName, ncol(common.matrix));
    cls.lbl <- dataSet$cls;
    
    for(i in 2:length(sel.nms)){
      dataName <- sel.nms[i];
      dataSet <- readRDS(dataName);
      ndat <- dataSet$data[shared.nms, ];
      nms.vec = c(nms.vec, rownames(dataSet$data));
      smps.vec = c(smps.vec, colnames(dataSet$data));
      plot.ndat <- t(scale(t(ndat)));
      common.matrix <- cbind(common.matrix, ndat);
      data.lbl <- c(data.lbl, rep(dataName, ncol(dataSet$data[,])));
      cls.lbl <- c(cls.lbl, dataSet$cls);
    }
    cls.lbl <- factor(cls.lbl);
    levels(cls.lbl) <- lvls;
    
    smps.nms = colnames(common.matrix)
    if(length(unique(smps.nms)) != length(smps.nms)){
      data.nb = length(unique(data.lbl));
      data.vec = vector()
      for(i in 1:data.nb){
        data.vec[i] = paste0("d", i);
      }
      levels(data.lbl) = data.vec;
      colnames(common.matrix) = make.unique(paste(data.vec, smps.nms, sep="_"));
      
      dataSet <- readRDS(sel.nms[1]);
      colnames(dataSet$data) = paste("d1", colnames(dataSet$data), sep="_");
      RegisterData(dataSet);
      
      for(i in 2:length(sel.nms)){
        dataSet <- readRDS(sel.nms[i]);
        colnames(dataSet$data) = paste0("d",i,"_",colnames(dataSet$data));
        # check if class label is consistent
        RegisterData(dataSet);
      }
      smps.vec <<- smps.nms;
      current.msg <<- paste("Duplicated sample names detected, samples have been renamed to make them unique.");
    }
    
    # note: index by entrez, gene symbol DO have duplicates
    rownames(common.matrix) <- shared.nms;
    
    # resort data, first on data.lbl, then on class lbl
    ord.inx <- order(data.lbl, cls.lbl);
    common.matrix <- data.matrix(common.matrix[,ord.inx]);
    cls.lbl <- cls.lbl[ord.inx];
    data.lbl <- data.lbl[ord.inx];
    smps.vec <- smps.vec[ord.inx];
    nms.vec = unique(nms.vec)
    nms.vec <<- nms.vec
    smps.vec <<- smps.vec
  }
  
  if(ncol(common.matrix) > 1000){  # max sample number allow 1000
    current.msg <<- paste("Total combined sample #:", ncol(common.matrix), "(exceed the limit: 1000!)");
    return(0);
  }
  
  # save the meta-dataset
  res <- data.frame(colnames(common.matrix), cls.lbl, data.lbl, t(common.matrix));
  colnames(res) <- c('#NAME', '#CLASS.condition', paste('#CLASS.dataset',id.type, data.org, sep="."), rownames(common.matrix));
  write.table(t(res), sep="\t", file="NetworkAnalyst_merged_data.txt", col.names=F, quote=FALSE);
  
  # need to set up the data for plotting (boxplot, heatmpa) so 
  # we need to scale row for each dataset in order to elimiate the maganitude difference 
  plot.matrix <- matrix(NA, nrow=nrow(common.matrix), ncol=ncol(common.matrix));
  rownames(plot.matrix) <- rownames(common.matrix);
  colnames(plot.matrix) <- colnames(common.matrix);
  for(i in 1:length(sel.nms)){
    data.inx <- data.lbl == sel.nms[i];
    plot.matrix[,data.inx] <- t(scale(t(common.matrix[,data.inx])));
  }
  
  # if entrez, get symbols for display
  shared.nms <- rownames(common.matrix);
  if(id.type == "entrez"){ 
    symbols <- doEntrez2SymbolMapping(shared.nms);
  }else{ # display itself
    symbols <- shared.nms;
  }
  names(symbols) <- shared.nms;
  
  inmex.meta <- list(data=common.matrix,
                     plot.data=plot.matrix,
                     id.type = id.type,
                     gene.symbls = symbols,
                     cls.lbl=factor(cls.lbl),
                     data.lbl=data.lbl);
  
  saveRDS(inmex.meta, "inmex_meta.rds");
  smps.vec <<- colnames(common.matrix);
  
  # setup common stats gene number, smpl number, grp info
  current.msg <<- paste("Sample #:", ncol(inmex.meta$data),
                        "Common ID #:", nrow(inmex.meta$data), 
                        "Condition:", paste(levels(inmex.meta$cls.lbl), collapse=" vs. "));
  return(1);
}

PerformMetaDeAnal <- function(){ 
  inmex.meta <- readRDS("inmex_meta.rds");   
  data.lbl = inmex.meta$data.lbl
  allmat = matrix("NA", nrow=length(nms.vec), ncol=length(data.lbl))
  rownames(allmat) = nms.vec
  colnames(allmat) = colnames(inmex.meta$data)
  
  include.inx <- mdata.all==1;
  if(sum(include.inx) == 0){
    current.msg <<-"No dataset is selected for analysis!";
    print(current.msg);
    return(0);
  }
  
  if(meta.upload){
    dat = inmex.meta$data 
    datasets = levels(inmex.meta$data.lbl)
    for(i in 1:length(datasets)){
      dt = datasets[i];
      inx = which(inmex.meta$data.lbl == dt)
      ind_dat = inmex.meta$data[,inx]
      cols <- colnames(allmat)[colnames(allmat) %in% colnames(ind_dat)]
      rows <- rownames(allmat)[rownames(allmat) %in% rownames(ind_dat)]
      norm.dat = t(apply(ind_dat, 1, function(x){as.numeric(cut(x, breaks=30))}))
      norm.dat = as.matrix(norm.dat)
      colnames(norm.dat) = colnames(ind_dat)
      allmat[rows, cols] <- norm.dat[rows, cols]
    }
  }else{
    sel.nms <- names(mdata.all)[include.inx];
    for(i in 1:length(sel.nms)){
      dataName <- sel.nms[i];
      dataSet <- readRDS(dataName);
      dataSet$data.orig <- NULL;
      cols <- colnames(allmat)[colnames(allmat) %in% colnames(dataSet$data)]
      rows <- rownames(allmat)[rownames(allmat) %in% rownames(dataSet$data)]
      norm.dat = t(apply(dataSet$data, 1, function(x){as.numeric(cut(x, breaks=30))}))
      norm.dat = as.matrix(norm.dat)
      colnames(norm.dat) = colnames(dataSet$data)
      allmat[rows, cols] <- norm.dat[rows, cols]
    }
  }
  
  saveRDS(allmat, "allmat.rds");
  performedDE <<- TRUE;
  PerformEachDEAnal(meta.upload);
}


# perform DE analysis on individual data (w.r.t common matrix)
# to be used/compared in the later analysis, with p-val Inf so that
# de can be adjusted based on user specified in meta later
PerformEachDEAnal <- function(is.meta=F){
  inmex.ind <- list();
  inmex.meta <- readRDS("inmex_meta.rds");
  sel.nms <- names(mdata.all)[mdata.all==1];
  if(is.meta){
    for(i in 1:length(sel.nms)){
      dataSet= list()
      dataName <- sel.nms[i];
      sel.inx <- inmex.meta$data.lbl == dataName;
      group <- factor(inmex.meta$cls.lbl[sel.inx]); # note regenerate factor to drop levels 
      data <- inmex.meta$data[, sel.inx];
      
      # update data set
      dataSet$type <- "array";
      dataSet$name <- dataName;
      group <- factor(inmex.meta$cls.lbl[sel.inx]); # note regenerate factor to drop levels 
      dataSet$cls <- group;
      res.limma <- PerformLimma(data, group);
      
      # save the limma fit object for meta-analysis (such as "dataSet1.fit.obj")
      saveRDS(res.limma$fit.obj, file=paste(dataName, "fit.obj", sep="."));
      
      res.all <- GetLimmaResTable(res.limma$fit.obj);
      saveRDS(res.all, "meta.resTable.rds");
      res.mat <- cbind(logFC=res.all$logFC, Pval = res.all$adj.P.Val);
      
      rownames(res.mat) <- rownames(res.all);
      inmex.ind[[dataName]] <- res.mat;
      
      #register sig one
      sig.inx <- res.mat[,2]<=GlobalCutOff$BHth;
      #sig.inx <- res.mat[,2]<=dataSet$pval;
      dataSet$sig.mat <- res.mat[sig.inx,];
      RegisterData(dataSet);
      
      # clean up
      rm(dataSet, res.all);
      gc();
    }
  }else{
    for(i in 1:length(sel.nms)){
      dataName <- sel.nms[i];
      sel.inx <- inmex.meta$data.lbl == dataName;
      group <- factor(inmex.meta$cls.lbl[sel.inx]); # note regenerate factor to drop levels 
      data <- inmex.meta$data[, sel.inx];
      
      dataSet <- readRDS(dataName);
      grp.lvl <- levels(dataSet$cls);
      
      # update data set
      dataSet$type <- "array";
      group <- factor(inmex.meta$cls.lbl[sel.inx], levels=grp.lvl, ordered=T); # note regenerate factor to drop levels 
      dataSet$cls <- group;
      
      res.limma <- PerformLimma(data, group);
      
      # save dataSet object for meta-analysis
      saveRDS(res.limma$fit.obj, file=paste(dataName, "fit.obj", sep="."));
      
      res.all <- GetLimmaResTable(res.limma$fit.obj);
      saveRDS(res.all, "meta.resTable.rds");
      
      res.mat <- cbind(logFC=res.all$logFC, Pval = res.all$adj.P.Val);
      
      rownames(res.mat) <- rownames(res.all);
      inmex.ind[[dataName]] <- res.mat;
      
      #sig.inx <- res.mat[,2]<=GlobalCutOff$BHth;
      sig.inx <- res.mat[,2]<=dataSet$pval;
      dataSet$sig.mat <- res.mat[sig.inx,];
      RegisterData(dataSet);
      # clean up
      rm(dataSet, res.all);
      gc();
    }
  }
  inmex.ind <<- inmex.ind;
}

# Meta-analysis combining effect size
PerformMetaEffectSize<- function(method="rem", BHth=0.05){
  if(!performedDE){
    PerformMetaDeAnal();
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  inmex.method <<- "effectsize";
  meta.mat <<- meta.stat <<- NULL;
  
  sel.nms <- names(mdata.all)[mdata.all==1];
  nbstudies <- length(sel.nms);
  listgd<-vector("list", (nbstudies+3));
  ES<-array(dim=c(nrow(inmex.meta$data),4,nbstudies));
  cls.lvls <- levels(as.factor(inmex.meta$cls.lbl));
  
  for (i in 1:nbstudies){
    data.nm <- sel.nms[i];
    dataSet <- readRDS(data.nm);
    
    fit.obj.nm <- paste(data.nm, "fit.obj", sep=".");
    fit2i <- readRDS(fit.obj.nm);
    
    pvals <- p.adjust(fit2i$p.value,method="BH");
    listgd[[i]]=which(pvals<=BHth);
    
    n1i=length(which(dataSet$cls==cls.lvls[1]));
    n2i=length(which(dataSet$cls==cls.lvls[2]));
    ES[,,i]=effectsize(fit2i$t,((n1i*n2i)/(n1i+n2i)),(fit2i$df.prior+fit2i$df.residual));
  }
  
  #only unbiased; for biased effect sizes, add ES[,1,],ES[,2,]
  listgd[[(nbstudies+1)]]=unique(unlist(listgd[1:nbstudies]))
  restempdirect=combineES(ES[,3,],ES[,4,],BHth, method);
  
  pooled.ef <- restempdirect[[3]];
  wt.mat<- restempdirect[[4]]; # one column for each study
  
  listgd[[(nbstudies+2)]]=restempdirect$DEindices
  listgd[[(nbstudies+3)]]=restempdirect$TestStatistic
  names(listgd)=c(paste("study",1:nbstudies,sep=""),"AllIndStudies","Meta","TestStatistic")  
  
  es.mat <- matrix(0, nrow=nrow(inmex.meta$data), ncol=2);
  es.mat[,1] <- pooled.ef;
  es.mat[,2] <- p.adjust(2*(1-pnorm(abs(listgd[[(nbstudies+3)]]))), method="BH");
  
  colnames(es.mat) <- c("CombinedES","Pval");
  rownames(es.mat) <- rownames(inmex.meta$data);
  #allMeta.mat <<- es.mat;
  saveRDS(es.mat, "allMeta.mat.rds");
  
  # now keep only genes with at least on sig (in one study or meta analysis)
  inx <- union(listgd[[(nbstudies+1)]], listgd[[(nbstudies+2)]]);
  es.mat <- es.mat[inx,];
  
  #sort
  ord.inx <- order(abs(es.mat[, "Pval"]), decreasing = F);
  es.mat<-signif(es.mat[ord.inx,],5);
  
  sig.inx <- which(es.mat[,"Pval"]<=BHth);
  meta.mat <<- es.mat[sig.inx, ];
  meta.mat.all <<- es.mat;
  SetupMetaStats(BHth);
  return(length(sig.inx));
}

# combining p values based on Fisher's or Stouffer
PerformPvalCombination <- function(method="stouffer", BHth=0.05){
  if(!performedDE){
    PerformMetaDeAnal();
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  inmex.method <<- "metap";
  meta.mat <<- meta.stat <<- NULL;
  sel.nms <- names(mdata.all)[mdata.all==1];
  
  classes <- list();
  nbstudies <- length(sel.nms);
  listgd=vector("list", (nbstudies+3));
  
  for (i in 1:nbstudies){
    data.nm <- sel.nms[i];
    dataSet <- readRDS(data.nm);
    classes[[i]] <- dataSet$cls; 
    
    fit.obj.nm <- paste(data.nm, "fit.obj", sep=".");
    fit2i <- readRDS(fit.obj.nm);
    
    pvals <- p.adjust(fit2i$p.value,method="BH");
    listgd[[i]]=which(pvals<=BHth);
    
    #recalculate moderated one sided p
    p1sidedLimma=pt(fit2i$t,df=(fit2i$df.prior+fit2i$df.residual))
    assign(paste("p1sidedLimma",i,sep=""), p1sidedLimma)
  }
  
  names(classes) <- sel.nms;
  tempvec=paste("p1sidedLimma",1:nbstudies,sep="");
  
  lsinglep=lapply(tempvec,FUN=function(x) get(x,inherits=TRUE));
  nrep=unlist(lapply(classes,FUN=function(x)length(x)));
  listgd[[(nbstudies+1)]]=unique(unlist(listgd[1:nbstudies]));
  
  restempdirect=combinePvals(lsinglep,nrep,BHth,method);
  
  listgd[[(nbstudies+2)]]=restempdirect$DEindices
  listgd[[(nbstudies+3)]]=restempdirect$CombinedP
  names(listgd)=c(paste("study",1:nbstudies,sep=""),"AllIndStudies","Meta","CombinedP");  
  
  pc.mat <- cbind(CombinedTstat=restempdirect$CombinedStat, CombinedPval=restempdirect$CombinedP);
  rownames(pc.mat) <- rownames(inmex.meta$data);
  saveRDS(pc.mat, "allMeta.mat.rds");
  
  # now keep only genes with at least on sig (in one study or meta analysis)
  inx <- union(listgd[[(nbstudies+1)]], listgd[[(nbstudies+2)]]);
  pc.mat <- pc.mat[inx,];
  
  #sort
  ord.inx <- order(pc.mat[, "CombinedPval"], decreasing = F);
  pc.mat<-signif(pc.mat[ord.inx,],5);
  
  sig.inx <- which(pc.mat[, "CombinedPval"]<=BHth);
  meta.mat <<- pc.mat[sig.inx, ];
  meta.mat.all <<- pc.mat;
  SetupMetaStats(BHth);
  return(length(sig.inx));
}

# diff used for direction, not selection
PerformVoteCounting <- function(BHth = 0.05, minVote){
  if(!performedDE){
    PerformMetaDeAnal();
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  inmex.method <<- "votecount";
  DE.vec <<- NULL; # store entrez id from meta-analysis for GO
  meta.mat <<- meta.stat <<- NULL;
  sel.nms <- names(mdata.all)[mdata.all==1];
  # first create a matrix to stall the result
  # row for each feature and col for each dataset uploaded
  vc.mat <- matrix(0, nrow=nrow(inmex.meta$data), ncol=length(sel.nms)+1);
  shared.ids <- rownames(inmex.meta$data);
  for(i in 1:length(inmex.ind)){
    res.mat <- inmex.ind[[i]];
    res.mat <- res.mat[shared.ids, ];
    
    #note in meta-analysis should consider directions
    # use logFC for this purpose 
    # consider upregulated
    hit.up.inx <- res.mat[,1]> 0 & res.mat[,2] <= BHth;
    up.vote <- as.numeric(hit.up.inx);
    
    # consider downregulated
    hit.dn.inx <- res.mat[,1] < 0 & res.mat[,2] <= BHth;
    dn.vote <- -as.numeric(hit.dn.inx);
    
    vc.mat[,i] <- up.vote + dn.vote;
  }
  
  # total score (votes for each direction)
  vc.mat[,length(sel.nms)+1] <- apply(vc.mat, 1, sum);
  colnames(vc.mat) <- c(paste("Vote", substring(sel.nms,0, nchar(sel.nms)-4)), "VoteCounts");
  rownames(vc.mat) <- rownames(inmex.meta$data);
  
  # compute at least one vote (no direction)
  vote.any <- apply(abs(vc.mat), 1, sum)
  vote.any.inx <- vote.any > 0;
  
  # return results with at least one vote
  vc.mat <- vc.mat[vote.any.inx, ];
  
  #sort
  ord.inx <- order(abs(vc.mat[, "VoteCounts"]), decreasing = T);
  vc.mat <- vc.mat[ord.inx, "VoteCounts", drop=F];
  
  sig.inx <- abs(vc.mat[,"VoteCounts"]) >= minVote;
  meta.mat <<- vc.mat;
  meta.mat.all <<- vc.mat;
  SetupMetaStats(BHth);
  return(length(sig.inx));
}


# This approach directly merge all data sets
# and analyze it as a single data
PerformMetaMerge<-function(BHth=0.05){
  inmex.method <<- "merge";
  meta.mat <<- meta.stat <<- NULL;
  inmex.meta <- readRDS("inmex_meta.rds");
  # prepare for meta-stats
  # calculate sig genes for individual analysis
  shared.names <- rownames(inmex.meta$data);
  
  res.limma <- PerformLimma(inmex.meta$data, as.factor(inmex.meta$cls.lbl));
  res.all <- GetLimmaResTable(res.limma$fit.obj);
  
  ord.inx <- order(res.all$adj.P.Val, decreasing=F);
  dm.mat <- as.matrix(res.all[ord.inx,c("logFC", "adj.P.Val")]);
  colnames(dm.mat) <- c("CombinedLogFC", "Pval");
  
  sig.inx <- which(dm.mat[,"Pval"] <= BHth);
  meta.mat <<- dm.mat[sig.inx,];
  meta.mat.all <<- dm.mat
  SetupMetaStats(BHth);
  return(length(sig.inx));
}

GetMetaGeneIDType<-function(){
  inmex.meta <- readRDS("inmex_meta.rds");
  return(inmex.meta$id.type);
}

GetMetaResultGeneIDs<-function(){
  rnms <- rownames(as.matrix(meta.mat));# already sorted based on meta-p values
  if(length(rnms) > 500){
    rnms <- rnms[1:500];
  }
  return(rnms);
}

# note, due to limitation of get/post
# maximum gene symb for list is top 500
GetMetaResultGeneSymbols<-function(){
  ids <- rownames(as.matrix(meta.mat));
  if(length(ids) > 500){
    ids <- ids[1:500];
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  if(inmex.meta$id.type == "entrez"){ # row name gene symbols
    ids <- inmex.meta$gene.symbls[ids];
  }
  return(ids);
}

GetMetaResultPathSymbols<-function(){
  return(rownames(meta.mat.all));
}

GetMetaResultPathLinks <- function(){
  symbs <- rownames(meta.mat);
  ids <- current.setids[symbs];
  # set up links to genbank
  annots <- paste("<a href='http://pantherdb.org/panther/category.do?categoryAcc=", ids,
                  "' target='_blank'>", symbs, "</a>", sep="");
  return(annots);
}


GetMetaResultGeneIDLinks <- function(){
  ids <- rownames(as.matrix(meta.mat));
  if(length(ids) > 500){
    ids <- ids[1:500];
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  symbs <- inmex.meta$gene.symbls[ids];
  # set up links to genbank
  annots <- paste("<a href='http://www.ncbi.nlm.nih.gov/gene?term=", ids,
                  "' target='_blank'>", symbs, "</a>", sep="");
  return(annots);
}

GetMetaResultColNames<-function(){
  sel.nms <- names(mdata.all)[mdata.all==1];
  c(substring(sel.nms, 0, nchar(sel.nms)-4), colnames(meta.mat));
}

# single.type return logFC or p value for individual data analysis
GetMetaResultMatrix<-function(single.type="fc"){
  if(single.type == "fc"){
    meta.mat <- cbind(fc.mat, meta.mat);
  }else{
    meta.mat <- cbind(pval.mat, meta.mat);
  }
  # display at most 500 genes
  if(nrow(meta.mat) > 500){
    meta.mat <- meta.mat[1:500,]; # already sorted based on meta-p values
  }
  meta.mat <-signif(as.matrix(meta.mat), 5);
  meta.mat;
}

GetMetaStat<-function(){
  return (meta.stat$stat);
}

GetMetaStatNames<-function(){
  return (names(meta.stat$stat));
}

combinePvals <- function(pvalonesided,nrep,BHth=0.05, method) {
  listres=vector("list",3);
  nbstudies=length(pvalonesided);
  nbreptot=sum(nrep);
  if (nbreptot <2) {
    stop("Error: the argument \"nrep\" must be a vector with at least two values higher than 1")
  } 
  
  weight=sqrt(nrep/nbreptot);
  fstatmeta=function(g){
    vecptime=unlist(lapply(pvalonesided, FUN = function(x) x[g]));
    vec = qnorm(1 - vecptime);
    stattestg = sum(weight[1:length(pvalonesided)] * vec[1:length(pvalonesided)], na.rm = TRUE);
    stattestg;
  }
  
  fishersum <- function(pvec){
    return(sum(-2*log(pvec)))
  }
  
  if(method=="stouffer"){
    statpvalc=-unlist(lapply(rep(1:length(as.vector(pvalonesided[[1]])), 1), function(x) fstatmeta(x)));
    rpvalpvalc=2*(1-pnorm(abs(statpvalc)));
  }else{ # fisher
    data <- data.frame(pvalonesided);
    #data[data == 0] <- 1*10^-10;
    
    #note, p value are calculated for one side
    # pt (lower.tail=T by default) which tests if group A < group B
    # for one side
    fsum1 <- apply(data, 1, fishersum);
    rpvalpvalc1 = 1-pchisq(fsum1, df=(ncol(data)*2));
    
    # for the other side
    data <- 1-data;
    fsum2 <- apply(data, 1, fishersum);
    rpvalpvalc2 = 1-pchisq(fsum2, df=(ncol(data)*2));
    
    # report the min of two direction calculation
    rpvalpvalc <- pmin(rpvalpvalc1, rpvalpvalc2);
    
    # adding direction information
    statpvalc <- pmax(fsum1, fsum2);
    # if A<B sig, then it should be negative 
    statpvalc[statpvalc == fsum1]<- -statpvalc[statpvalc == fsum1];
  }
  
  rpvalpvalc <- p.adjust(rpvalpvalc,method="BH");
  res=which(rpvalpvalc<=BHth);
  listres[[1]]=res
  listres[[2]]=statpvalc;
  listres[[3]]=rpvalpvalc
  names(listres)=c("DEindices", "CombinedStat", "CombinedP")
  listres
}

#combine effect size
combineES <- function (ES, varES, BHth = 0.05, method){
  if(method == "rem"){
    useREM = TRUE;
  }else{
    useREM = FALSE;
  }
  
  num.studies <- dim(ES)[2];
  
  Qvals <- f.Q.NA(ES, varES)
  if (useREM) {
    varES <- varES + tau2.NA(Qvals, num.studies, my.weights = 1/varES)
  }
  wt <- 1/varES
  MUvals <- rowSums(ES * wt, na.rm = TRUE)/rowSums(wt, na.rm = TRUE)
  MUsES <- sqrt(abs(1/rowSums(wt, na.rm = TRUE)))
  zSco <- MUvals/MUsES
  rpvalESc = 2 * (1 - pnorm(abs(zSco)))
  res = which(p.adjust(rpvalESc, method = "BH") <= BHth);
  listres <- list();
  listres[[1]] = res
  listres[[2]] = zSco
  listres[[3]] = MUvals; # pool effect size
  listres[[4]] = wt; # wt for each studies, it is matrix with one column for each studies
  names(listres) = c("DEindices", "TestStatistic", "PooledEffectSize", "Weights")
  listres
}

# first two biased ES, var, last two unbiased

f.Q.NA = function(dadj, varadj) {
  w <- 1/varadj
  tmp1 <- w * dadj
  mu <- rowSums(tmp1, na.rm = TRUE)/rowSums(w, na.rm = TRUE)
  Q <- rowSums(w * (dadj - mu)^2, na.rm = TRUE)
}

tau2.NA <- function(Q, num.studies, my.weights) {
  vwts <- rowSums(my.weights, na.rm = TRUE)
  tmp2 <- rowSums(my.weights^2, na.rm = TRUE)
  tau2 <- pmax(0, (Q - (num.studies - 1))/(vwts - tmp2/vwts))
  return(tau2)
}

# prepare data for heatmap plotting include
# 1. meta info (top bars)
# 2. expression matrix
# 3. function annotation (left bars)
# all matrix will be combined, 
# 1 and 2 separated by a row of 'null' 
# 3 and 1+2 separated by a column of 'null'
PrepareMetaHeatmapJSON <- function(){
  gene.vec <- rownames(meta.mat);
  datanm.vec <- names(mdata.all)[mdata.all==1];
  
  inmex.meta <- readRDS("inmex_meta.rds");
  dat.inx <- inmex.meta$data.lbl %in% datanm.vec;
  dat <- inmex.meta$plot.data[gene.vec, dat.inx, drop=F]; 
  
  # scale each gene for each dataset
  dat <- t(scale(t(dat)));
  
  # now need to remove na or constant rows
  dat <- na.omit(dat);
  # check for columns with all constant (var =0)
  varCol <- apply(dat, 1, var, na.rm=T);
  constCol <- (varCol == 0 | is.na(varCol));
  dat <- dat[!constCol, ];
  
  anot.res <- list();
  ids <- rownames(dat);
  if(inmex.meta$id.type == "entrez"){   
    anot.res <- doEntrezIDAnot(ids);
  }else{ # no annotation, then use the default feature ID
    anot.res$gene_id <- anot.res$symbol <- anot.res$name <- ids; 
  }
  
  data.nms <- as.factor(inmex.meta$data.lbl[dat.inx]);
  cls.lbls <- as.factor(inmex.meta$cls.lbl[dat.inx]);
  
  # setup annotation info
  data.nms <- as.character(data.nms);
  datasets <- substr(as.character(data.nms), 0, nchar(data.nms)-4);
  
  sel.nms <- names(mdata.all)[mdata.all==1];
  if(length(sel.nms) > 1){
    annotation <- data.frame(class= cls.lbls, dataset = as.factor(datasets));
  }else{ # single data
    annotation <- data.frame(class= cls.lbls);
  }
  rownames(annotation) <- colnames(dat);
  
  ####
  sig.ids <- rownames(dat);
  if(inmex.method != "votecount"){
    stat.pvals <- as.numeric(meta.mat[,2]);
    stat.fc = as.numeric(meta.mat[,1]);
  }else{
    stat.pvals <- as.numeric(meta.mat[,1]);
    stat.fc = as.numeric(meta.mat[,1]);
  }
  
  orig.smpl.nms <- colnames(dat);
  orig.gene.nms <- rownames(dat);
  
  # do clustering and save cluster info
  if(nrow(dat)> 1){
    dat.dist <- dist(dat);
    gene.ward.ord <- hclust(dat.dist, "ward.D")$order;
    gene.ward.rk <- match(orig.gene.nms, orig.gene.nms[gene.ward.ord]);
    gene.ave.ord <- hclust(dat.dist, "ave")$order;
    gene.ave.rk <- match(orig.gene.nms, orig.gene.nms[gene.ave.ord]);
    gene.single.ord <- hclust(dat.dist, "single")$order;
    gene.single.rk <- match(orig.gene.nms, orig.gene.nms[gene.single.ord]);
    gene.complete.ord <- hclust(dat.dist, "complete")$order;
    gene.complete.rk <- match(orig.gene.nms, orig.gene.nms[gene.complete.ord]);
    
    dat.dist <- dist(t(dat));
    smpl.ward.ord <- hclust(dat.dist, "ward.D")$order;
    smpl.ward.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ward.ord])
    smpl.ave.ord <- hclust(dat.dist, "ave")$order;
    smpl.ave.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.ave.ord])
    smpl.single.ord <- hclust(dat.dist, "single")$order;
    smpl.single.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.single.ord])
    smpl.complete.ord <- hclust(dat.dist, "complete")$order;
    smpl.complete.rk <- match(orig.smpl.nms, orig.smpl.nms[smpl.complete.ord])
  }else{
    # force not to be single element vector which will be scaler
    stat.pvals <- matrix(stat.pvals);
    gene.ward.rk <- gene.ave.rk <- gene.single.rk <- gene.complete.rk <- matrix(1);
    smpl.ward.rk <- smpl.ave.rk <- smpl.single.rk <- smpl.complete.rk <- 1:ncol(dat);
  }
  
  gene.cluster <- list(
    pval = stat.pvals,
    fc = stat.fc, 
    ward = gene.ward.rk,
    average = gene.ave.rk,
    single = gene.single.rk,
    complete = gene.complete.rk
  );
  
  sample.cluster <- list(
    ward = smpl.ward.rk,
    average = smpl.ave.rk,
    single = smpl.single.rk,
    complete = smpl.complete.rk
  );
  
  # prepare meta info    
  # 1) convert meta.data info numbers
  # 2) match number to string (factor level)
  meta <- annotation;
  grps <- colnames(meta);
  nmeta <- meta.vec <- NULL;
  uniq.num <- 0;
  for (i in 1:ncol(meta)){
    cls <- meta[,i];
    grp.nm <- grps[i];
    meta.vec <- c(meta.vec, as.character(cls))
    # make sure each label are unqiue across multiple meta data
    ncls <- paste(grp.nm, as.numeric(cls)); # note, here to retain ordered factor
    nmeta <- c(nmeta, ncls);
  }
  
  # convert back to numeric 
  nmeta <- as.numeric(as.factor(nmeta))+99;
  unik.inx <- !duplicated(nmeta)   
  
  # get corresponding names
  meta_anot <- meta.vec[unik.inx]; 
  names(meta_anot) <- nmeta[unik.inx]; # name annotatation by their numbers
  
  nmeta <- matrix(nmeta, ncol=ncol(meta), byrow=F);
  colnames(nmeta) <- grps;
  
  # for each gene/row, first normalize and then tranform real values to 30 breaks 
  res <- t(apply(dat, 1, function(x){as.numeric(cut(x, breaks=30))}));
  
  # note, use {} will lose order; use [[],[]] to retain the order
  # single element vector will be converted to scalar, not array, need to prevent that
  gene.id = anot.res$symbol; if(length(gene.id) ==1) { gene.id <- matrix(gene.id) };
  gene.entrez = anot.res$gene_id; if(length(gene.entrez) ==1) { gene.entrez <- matrix(gene.entrez) };        
  gene.name = anot.res$name; if(length(gene.name) ==1) { gene.name <- matrix(gene.name) };
  
  json.res <- list(
    data.type = "array",
    gene.id = as.character(anot.res$symbol),
    gene.entrez = gene.entrez,
    gene.name = anot.res$name,
    gene.cluster = gene.cluster,
    sample.cluster = sample.cluster,
    sample.names = orig.smpl.nms,
    meta = data.frame(nmeta),
    meta.anot = meta_anot,
    data.lbl = inmex.meta$data.lbl,
    data = res
  );
  
  return(json.res);
}


##################################
# functions for estimating Cochran’s Q
##################################


#computes Cochran’s Q gene by gene
#dadj and varadj must be matrices, in which every study is a column,
#every row a gene
f.Q <- function(dadj, varadj){
  w<-1/varadj
  tmp1<-w*dadj
  mu<-rowSums(tmp1)/rowSums(w)
  Q<-rowSums(w*(dadj - mu)^2)
}

qc.metaDensity<- function(imgNm, dpi=72, format, factor){
  library("ggplot2")
  inmex.meta <- readRDS("inmex_meta.rds");
  dat = inmex.meta$data
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  dpi = as.numeric(dpi)
  
  df = data.frame(inmex.meta$data, stringsAsFactors = FALSE)
  df = stack(df)
  if(factor != "NA"){
    #Factor = as.vector(dataSet$meta.info[,factor])
  }else{
    #Factor = dataSet$meta.info[,1];
  }
  Factor=inmex.meta$data.lbl
  
  conv = data.frame(ind=colnames(inmex.meta$data), class=Factor)
  conv$ind=gsub("-", ".", conv$ind)
  df1 = merge(df, conv, by="ind")
  Cairo(file=imgNm, width=10, height=6, type=format, bg="white", dpi=dpi, unit="in");
  g =ggplot(df1, aes(x=values)) + geom_line(aes(color=class, group=ind), stat="density", alpha=0.3) + geom_line(aes(color=class), stat="density", alpha=0.6, size=1.5)
  print(g)
  dev.off();
}

effectsize <- function(tstat,ntilde,m){
  cm=gamma(m/2)/(sqrt(m/2)*gamma((m-1)/2))
  d=tstat/sqrt(ntilde)
  dprime=cm*d
  terme1=m/((m-2)*ntilde)
  vard=terme1+d^2*(terme1*ntilde-1/cm^2)
  vardprime=cm^2*(terme1+dprime^2*(terme1*ntilde-1/cm^2))
  result=cbind(d,vard,dprime,vardprime)
  colnames(result)=c("d","vard","dprime","vardprime")
  result
}

#For internal GSEA of meta-analysis of gene sets

PerformMetaPathCombine <- function(name, netNm, method, lib, mType, BHth=0.05){
  if(!performedDE){
    PerformMetaDeAnal();
  }
  library(fgsea)
  curr.geneset <- LoadEnrLib(lib);
  inmex.method <<- "effectsize";
  meta.stat = "null";
  
  allMeta.mat <- readRDS("allMeta.mat.rds");
  #allMeta.mat[allMeta.mat[,2]==0] <- 1e-20
  rankedVec = as.vector(allMeta.mat[,1])*sign(allMeta.mat[,1]);
  
  names(rankedVec) = rownames(allMeta.mat);
  rankedVec = sort(rankedVec)
  rankedVec = rankedVec[unique(names(rankedVec))]
  fgseaRes <- fgsea(pathways = curr.geneset, 
                    stats = rankedVec,
                    minSize=1,
                    maxSize=1000,
                    nperm=5000)
  fgseaRes = fgseaRes[!duplicated(fgseaRes$pathway),]
  fgseaRes = fgseaRes[,c("size","ES", "NES","padj", "pathway", "pval")]
  fgseaRes=fgseaRes[order(-abs(fgseaRes$ES)),]
  fgseaRes=fgseaRes[order(fgseaRes$pval),] 
  
  fgseaRe = data.frame(fgseaRes)
  rownames(fgseaRes) = fgseaRes$pathway
  es.mat = as.matrix(fgseaRes[,c("ES","padj")]);
  colnames(es.mat) = c("EnrichmentScore","Pvalue")
  rownames(es.mat) = fgseaRes$pathway
  sig.inx <- which(es.mat[, "Pvalue"]<=BHth);
  if(length(sig.inx)<10){
    sig.inx = c(1:10)
  }
  metaset.mat <<- es.mat[sig.inx,];
  metaset.mat.all <<- fgseaRes[sig.inx,]
  ii = SetupMetaGSEAStats(name, netNm, BHth, mType, curr.geneset,lib);
  if(ii == 1){
    return(length(sig.inx));
  }else{
    return(0);
  }
}


PlotGShmMeta <-function(cmpdNm, IDs){
  ids = unlist(strsplit(IDs, "; "));
  cmpdNm <- gsub(" ", "_",  cmpdNm);
  cmpdNm <- gsub("/", "_",  cmpdNm);
  inmex.meta <- readRDS("inmex_meta.rds");
  subset = dataSet$data.norm[which(doEntrez2SymbolMapping(rownames(inmex.meta$data)) %in% ids),]
  if(length(subset)<1){
    subset = dataSet$data.norm[which(rownames(inmex.meta$data) %in% ids),]
  }
  json.res <- list(
    data=subset,
    ids=doEntrez2SymbolMapping(rownames(subset)),
    entrez = rownames(subset)
  )
  library(RJSONIO);
  json.mat <- RJSONIO::toJSON(json.res, .na='null', pretty=TRUE);
  json.nm <- paste(cmpdNm,"_hm", ".json", sep="");
  
  sink(json.nm);
  cat(json.mat);
  sink();
  
  return(json.nm)
}

PlotMetaHm <-function(cmpdNm){
  
  allmat = readRDS("allmat.rds")
  current.geneset <- readRDS("current_geneset.rds")
  ids <- current.geneset[[cmpdNm]];
  subset = allmat[which(rownames(allmat) %in% ids),]
  if(length(subset)<1){
    subset = allmat[which(rownames(allmat) %in% ids),]
  }
  
  if(inmex.method %in% c("effectsize", "merge")){
    dims = dim(subset)
    rnms = rownames(subset)
    cnms = colnames(subset) 
    m <- mapply(subset, FUN=as.numeric)
    subset <- matrix(data=m, ncol=dims[2], nrow=dims[1])
    rownames(subset) = rnms
    colnames(subset) = cnms
    subset = subset[complete.cases(subset), ];
  }
  
  json.res <- list(
    data=subset,
    ids=doEntrez2SymbolMapping(rownames(subset)),
    entrez = rownames(subset)
  )
  library(RJSONIO);
  json.mat <- RJSONIO::toJSON(json.res, .na='null', pretty=TRUE);
  json.nm <- paste(cmpdNm,"_hm", ".json", sep="");
  
  sink(json.nm);
  cat(json.mat);
  sink();
  
  return(json.nm)
}

PlotMetaPhm <-function(cmpdNm, dpi=72){
  allmat = readRDS("allmat.rds");
  
  fileNm = paste("Path_", cmpdNm, ".png", sep="");
  current.geneset <- readRDS("current_geneset.rds")
  ids=current.geneset[[cmpdNm]];
  subset = allmat[which(rownames(allmat) %in% ids),]
  if(length(subset)<1){
    subset = allmat[which(rownames(allmat) %in% ids),]
  }
  library(RColorBrewer);
  library(pheatmap)
  dims = dim(subset)
  rnms = rownames(subset)
  cnms = colnames(subset) 
  m <- mapply(subset, FUN=as.numeric)
  subset <- matrix(data=m, ncol=dims[2], nrow=dims[1])
  rownames(subset) = rnms
  colnames(subset) = cnms
  if(inmex.method %in% c("effectsize", "merge")){
    subset = subset[complete.cases(subset), ];
  }
  
  subset[is.na(subset)] <- 0
  my_palette <- colorRampPalette(c("green", "black", "red"))(n = 30)
  my_palette = c("#d3d3d3", my_palette)
  bk2 = unique(c(seq(0,0.9999, length=2), 1, seq(2,30, length=29)));
  inmex.meta <- readRDS("inmex_meta.rds");
  ann = data.frame(Class=inmex.meta$cls.lbl, Dataset=inmex.meta$data.lbl)
  rownames(ann) = colnames(subset)
  Cairo(file=fileNm, width=800, height=700, type="png", bg="white",unit="px",dpi=72)
  hm = pheatmap(subset, color = my_palette, breaks = bk2, annotation_col=ann, show_rownames = FALSE, main = cmpdNm, border_color=NA)
  print(hm)
  dev.off()
}

# for gene set-level meta-analysis

SetupMetaGSEAStats <- function(name, netNm, BHth, mType, curr.geneset, lib){
  
  inmex.de <- list();
  allmat <- readRDS("allMeta.mat.rds");
  allmat.vec <- rownames(allmat);
  meta.mat = metaset.mat
  metade.genes <- rownames(meta.mat);
  pval.mat <- fc.mat <- matrix(nrow=nrow(meta.mat), ncol=sum(mdata.all==1));
  fc.list <- split(rep(" ", length(allmat.vec)), allmat.vec);
  
  
  current.geneset = curr.geneset[!duplicated(names(curr.geneset))]
  inx =names(current.geneset) %in% rownames(meta.mat)  ;
  
  resTable = meta.mat
  current.mset = current.geneset[inx];
  
  inmex.meta <- readRDS("inmex_meta.rds");
  
  ora.vec <- rownames(inmex.meta$data)
  ora.nms <- doEntrez2SymbolMapping(rownames(inmex.meta$data))
  
  hits.query <- lapply(current.mset, 
                       function(x) {
                         unique(ora.nms[ora.vec%in%unlist(x)]);
                       });  
  saveRDS(hits.query, "hits_query.rds");
  
  set.num = unlist(lapply(current.mset, function(x){length(unique(x))}), use.names=TRUE);
  names(hits.query) <- names(current.mset);
  hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=TRUE);
  
  vote.bool = "false"
  meta.matcolinx = 2;
  enr.score = "NA"   
  
  padj <- p.adjust(as.vector(meta.mat[,meta.matcolinx]),method="BH");
  if(mType == "network"){
    json.res <- list(
      fun.anot = hits.query,
      fun.ids = as.vector(rownames(meta.mat)),
      fun.pval = as.vector(meta.mat[,meta.matcolinx]),
      fun.padj = padj,
      hit.num = hit.num,
      total= set.num
    );
  }else{
    json.res <- list(
      hits = hit.num,
      total= set.num,
      enr.pval= as.vector(meta.mat[,meta.matcolinx]),
      enr.padj= padj,
      enr.names= as.vector(rownames(meta.mat)),
      cls.lbl=inmex.meta$cls.lbl,
      smps.lbl=smps.vec,
      data.lbl = inmex.meta$data.lbl,
      path.lbl = rownames(meta.mat),
      enr.score = as.vector(meta.mat[,1]),
      isVote = vote.bool
    );
  }
  
  
  library(RJSONIO);
  json.mat <- RJSONIO::toJSON(json.res, .na='null', pretty=TRUE);
  json.nm <- paste0(name, ".json");
  
  sink(json.nm);
  cat(json.mat);
  sink();
  if(mType == "network"){
    res.mat<-matrix(0, nrow=length(names(current.mset)), ncol=4);
    rownames(res.mat)<-names(current.mset);
    colnames(res.mat)<-c("Total","Hits", "P.Value", "FDR");
    res.mat[,"Total"] = set.num
    res.mat[,"Hits"] = hit.num;
    res.mat[,"P.Value"] = meta.mat[,meta.matcolinx];
    res.mat[,"FDR"] = padj;
    enr.mat <<- res.mat
    res <- data.frame(Name=as.vector(rownames(meta.mat)), Total=set.num, Hits= hit.num, EnrichmentScore=as.vector(meta.mat[,1]), Pval=as.vector(meta.mat[,meta.matcolinx]), Padj = padj);
    list.genes <<- allmat.vec
    SetListNms();
    netnm <- paste0(netNm, ".json");
    PrepareEnrichNet(netNm, "meta", "mixed");
  }else{
    res.mat<-matrix(0, nrow=length(names(current.mset)), ncol=5);
    colnames(res.mat)<-c("Name", "Total","Hits", "P.Value", "FDR");
    res.mat[,"Name"] = names(current.mset);
    res.mat[,"Total"] = set.num
    res.mat[,"Hits"] = hit.num;
    res.mat[,"P.Value"] = meta.mat[,meta.matcolinx];
    res.mat[,"FDR"] = padj;
    write.csv(res.mat, file=paste("meta_sig_genesets_", lib, ".csv", sep=""), row.names=F);
  }
  return(1)
}

CalculateGsNet <- function(name, netNm, type, mType, db){
  res = PerformMetaPathCombine(name, netNm, "pval", db, mType,0.05)
  return(1)
}

##################################################
## R scripts for NetworkAnalyst 
## Various utility methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################
# given a text from input area: one or two-col entries (geneID and logFC)
# parse out return the a matrix containing the logFc, with rows named by gene ID
# if no 2nd col (logFC), the value will be padded by 0s
.parseListInput <- function(geneIDs){
  
  spl <- unlist(strsplit(geneIDs, "\\//")[1]);
  spl <- spl[unlist(lapply(spl,function(x){!x %in% ""}))]
  spl <- lapply(spl,function(x){gsub("\\/", "",x)})
  numOfLists <<- length(spl)
  dataList = list();
  inxU <- 0;
  for (i in 1:length(spl)){
    lines <- unlist(strsplit(spl[[i]], "\r|\n|\r\n")[1]);
    # remove the beginning & trailing space 
    lines <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", lines, perl=TRUE);
    if(substring(lines[1],1,1)=="#"){
      lines <- lines[-1];
    }
    gene.lists <- strsplit(lines, "\\s+");
    gene.mat <- do.call(rbind, gene.lists);
    
    if(dim(gene.mat)[2] == 1){ # add 0
      gene.mat <- cbind(gene.mat, rep(0, nrow(gene.mat)));
      current.msg <- "Only one column found in the list. Abundance will be all zeros. ";
    }else if(dim(gene.mat)[2] > 2){
      gene.mat <- gene.mat[,1:2];
      current.msg <- "More than two columns found in the list. Only first two columns will be used. ";
    }
    print(current.msg);
    
    rownames(gene.mat) <- gene.mat[,1];
    gene.mat <- gene.mat[,-1, drop=F];
    inxU <- inxU + 1;
    listInxU <<- paste0("datalist", inxU);
    gene.mat <- RemoveDuplicates(gene.mat, "mean", quiet=F); 
    good.inx <- !is.na(gene.mat[,1]);
    gene.mat <- gene.mat[good.inx, , drop=F];
    dataList[[i]] = gene.mat  
  }
  return(dataList)
}

GetSelListLength <- function(nm){
  if(dataSet$name != nm){
    dataSet = readRDS(nm);
  }
  return(length(dataSet$sig.mat));
}


GetColorSchema <- function(my.grps){
  # test if total group number is over 9
  my.grps = as.factor(my.grps);
  grp.num <- length(levels(my.grps));
  
  if(grp.num > 9){
    pal12 = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
              "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
              "#FFFF99", "#B15928");
    dist.cols <- colorRampPalette(pal12)(grp.num);
    lvs <- levels(my.grps);
    colors <- vector(mode="character", length=length(my.grps));
    for(i in 1:length(lvs)){
      colors[my.grps == lvs[i]] <- dist.cols[i];
    }
  }else{
    colors <- as.numeric(my.grps)+1;
  }
  return (colors);
}

GetExtendRange<-function(vec, unit=10){
  var.max <- max(vec);
  var.min <- min(vec);
  exts <- (var.max - var.min)/unit;
  c(var.min-exts, var.max+exts);
}

# given a data with duplicates, dups is the one with duplicates
RemoveDuplicates <- function(data, lvlOpt, quiet=T){
  
  all.nms <- rownames(data);
  colnms <- colnames(data);
  dup.inx <- duplicated(all.nms);
  dim.orig  <- dim(data);
  data <- apply(data, 2, as.numeric); # force to be all numeric
  dim(data) <- dim.orig; # keep dimension (will lost when only one item) 
  rownames(data) <- all.nms;
  colnames(data) <- colnms;
  if(sum(dup.inx) > 0){
    uniq.nms <- all.nms[!dup.inx];
    uniq.data <- data[!dup.inx,,drop=F];
    
    dup.nms <- all.nms[dup.inx];
    uniq.dupnms <- unique(dup.nms);
    uniq.duplen <- length(uniq.dupnms);
    
    for(i in 1:uniq.duplen){
      nm <- uniq.dupnms[i];
      hit.inx.all <- which(all.nms == nm);
      hit.inx.uniq <- which(uniq.nms == nm);
      
      # average the whole sub matrix 
      if(lvlOpt == "mean"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, mean, na.rm=T);
      }else if(lvlOpt == "median"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, median, na.rm=T);
      }else if(lvlOpt == "max"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, max, na.rm=T);
      }else{ # sum
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, sum, na.rm=T);
      }
    }
    if(!quiet){
      if(numOfLists == 1){
        current.msg <<- paste(current.msg, paste("A total of ", sum(dup.inx), " of duplicates were replaced by their ", lvlOpt, ".", sep=""), collapse="\n");
      }else{
        current.msg <<- paste(current.msg, paste0("<b>", listInxU, "</b> : ", length(data), " genes;"), collapse="\n");
      }
    }
    return(uniq.data);
  }else{
    if(!quiet){
      if(numOfLists == 1){
        current.msg <<- paste(current.msg, "All IDs are unique.", collapse="\n");
      }else{
        current.msg <<- paste(current.msg, paste0("<b>", listInxU, "</b> : ", length(data), " genes;"), collapse="\n");
      }
    }
    return(data);
  }
} 

# utils to remove from
# within, leading and trailing spaces
# remove /
ClearFactorStrings<-function(cls.nm, query){
  # remove leading and trailing space
  query<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", query, perl=TRUE);
  
  # kill multiple white space
  query <- gsub(" +","_",query);
  # remove non alphabets and non numbers 
  query <- gsub("[^[:alnum:] ]", "_", query);
  
  # test all numbers (i.e. Time points)
  chars <- substr(query, 0, 1);
  num.inx<- chars >= '0' & chars <= '9';
  if(all(num.inx)){
    query = as.numeric(query);
    nquery <- paste(cls.nm, query, sep="_");
    query <- factor(nquery, levels=paste(cls.nm, sort(unique(query)), sep="_"));
  }else{
    query[num.inx] <- paste(cls.nm, query[num.inx], sep="_");
    query <- factor(query);
  }
  return (query);
}

# borrowed from Hmisc
all.numeric <- function (x, what = c("test", "vector"), extras = c(".", "NA")){
  what <- match.arg(what)
  old <- options(warn = -1)
  on.exit(options(old));
  x <- sub("[[:space:]]+$", "", x);
  x <- sub("^[[:space:]]+", "", x);
  inx <- x %in% c("", extras);
  xs <- x[!inx];
  isnum <- !any(is.na(as.numeric(xs)))
  if (what == "test") 
    isnum
  else if (isnum) 
    as.numeric(x)
  else x
}

# utils to remove from
# within, leading and trailing spaces
# remove /
ClearStrings<-function(query){
  # remove leading and trailing space
  query<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", query, perl=TRUE);
  
  # kill multiple white space
  query <- gsub(" +",".",query);
  query <- gsub("/", ".", query);
  query <- gsub("-", ".", query);
  return (query);
}

# need to obtain the full path to convert (from imagemagik) for cropping images
GetBashFullPath<-function(){
  path <- system("which bash", intern=TRUE);
  if((length(path) == 0) && (typeof(path) == "character")){
    print("Could not find bash in the PATH!");
    return("NA");
  }
  return(path);
}

#######################################
### Utility Methods not for public call
########################################
# note, last two par only for STRING database
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

# private method for all sqlite queries
.query.sqlite <- function(db.con, statement, offline=TRUE){
  rs <- dbSendQuery(db.con, statement);
  res <- fetch(rs, n=-1); # get all records
  dbClearResult(rs);
  if(offline){
    dbDisconnect(db.con);
  }
  cleanMem();
  return(res);
}

cleanMem <- function(n=8) { for (i in 1:n) gc() }

###########
# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    capture.output(format(utils::object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  print(lapply(dataSet, object.size));
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

# shorthand
ShowMemoryUse <- function(..., n=30) {
  library(pryr);
  sink(); # make sure print to screen
  print(mem_used());
  print(sessionInfo());
  print(.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n));
  print(warnings());
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

LoadCurrentSet <- function() {
  currentSet <- readRDS("current_geneset.rds");
  return(currentSet);
}

PerformHeatmapEnrichment <- function(file.nm, fun.type, IDs){
  if(IDs=="NA"){
    if(anal.type=="onedata"){
      gene.vec <- rownames(dataSet$sig.mat);
    }else if(anal.type=="metadata"){
      gene.vec <- rownames(meta.mat);
    }else{
      gene.vec <- rownames(all.ent.mat);
    }
  }else{
    gene.vec <- unlist(strsplit(IDs, "; "));
  }
  sym.vec <- doEntrez2SymbolMapping(gene.vec);
  names(gene.vec) <- sym.vec;
  res <- PerformEnrichAnalysis(file.nm, fun.type, gene.vec);
  return(res);
}

PrepareEnrichNet<-function(netNm, type, overlapType){
  hits <-  enr.mat[,"Hits"];
  pvals <- enr.mat[,"P.Value"];
  require(igraph);
  require(reshape);
  pvalue <- pvals;
  id <- names(pvalue);
  current.geneset <- readRDS("current_geneset.rds");
  hits.query <- readRDS("hits_query.rds")
  hits.query <- hits.query[rownames(enr.mat)];
  geneSets <- hits.query;
  n <- length(pvalue);
  w <- matrix(NA, nrow=n, ncol=n);
  colnames(w) <- rownames(w) <- id;
  
  for (i in 1:n) {
    for (j in i:n) {
      w[i,j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]], overlapType)
    }
  }
  wd <- melt(w);
  wd <- wd[wd[,1] != wd[,2],];
  wd <- wd[!is.na(wd[,3]),];
  
  g <- graph.data.frame(wd[,-3], directed=F);
  if(type == "list"){
    g <- delete.edges(g, E(g)[wd[,3] < 0.3]);
  }else{
    g <- delete.edges(g, E(g)[wd[,3] < 0.3]);
  }
  idx <- unlist(sapply(V(g)$name, function(x) which(x == id)));
  
  # define local function
  my.normalize <- function(x){
    return((x- min(x)) /(max(x)-min(x)))
  }
  my.rescale <- function(x, from, to){
    (x - min(x)) / max(x - min(x)) * (to - from) + from
  }
  
  V(g)$color <- ComputeColorGradient(-log(my.normalize(pvalue) + min(pvalue/2)), "black", F, F);
  V(g)$colorw <- ComputeColorGradient(-log(my.normalize(pvalue) + min(pvalue/2)), "white", F, F);
  
  cnt <- hits;
  names(cnt) <- id;
  cnt2 <- cnt[V(g)$name];
  
  V(g)$size <- my.rescale(log(cnt2, base=10), 8, 32);
  
  # layout
  pos.xy <- layout.auto(g);
  
  # now create the json object
  nodes <- vector(mode="list");
  node.nms <- V(g)$name;
  node.sizes <- V(g)$size;
  node.cols <- V(g)$color;
  node.colsw <- V(g)$colorw;
  
  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id = node.nms[i],
      label=node.nms[i],
      size = node.sizes[i],
      true_size=node.sizes[i], 
      colorb=node.cols[i],
      colorw=node.colsw[i],
      posx = pos.xy[i,1],
      posy = pos.xy[i,2]
    );
  }
  
  edge.mat <- get.edgelist(g);
  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2]);
  
  # covert to json
  bedges = stack(hits.query);
  b.mat <- matrix(NA, nrow=nrow(bedges), ncol=2);
  b.mat[,1] = bedges[,"values"];
  b.mat[,2] = as.character(bedges[,"ind"]);
  b.mat = b.mat[complete.cases(b.mat),]
  colnames(b.mat) = c("source", "target");
  bg <- graph.data.frame(b.mat, directed=F);
  idx <- unlist(sapply(V(bg)$name, function(x) which(x == id)));
  cols <- color_scale("red", "#E5C494");
  
  V(bg)$color[V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(-log(pvalue), "black", F, F);
  V(bg)$colorw[V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(-log(pvalue), "white", F, F);
  node.nms <- V(bg)$name;
  if(anal.type == "onedata"){
    tbl = dataSet$resTable
    tbl = tbl[which(doEntrez2SymbolMapping(rownames(tbl)) %in% V(bg)$name),]
    expr.val = tbl[,selectedFactorInx];
    expvals = expr.val;
    names(expvals) = doEntrez2SymbolMapping(rownames(tbl))
    expvals = expvals[node.nms]
    V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T, F);
    V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T, F);
  }else if(anal.type == "genelist" && sum(all.prot.mat[,1]) != 0){
    tbl = all.prot.mat
    gene.nms = V(bg)$name[which(!V(bg)$name %in% rownames(enr.mat))]
    tbl = tbl[which(rownames(tbl) %in% gene.nms),]
    expr.val = tbl;
    expvals = expr.val
    expvals = expvals[node.nms]
    V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T, F);
    V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T, F);
    
  }else if(anal.type =="metadata"){
    tbl = meta.mat.all
    tbl = as.matrix(tbl[which(doEntrez2SymbolMapping(rownames(tbl)) %in% V(bg)$name),])
    expr.val = tbl[,1];
    expvals = expr.val;
    names(expvals) = doEntrez2SymbolMapping(rownames(tbl))
    expvals = expvals[node.nms]
    V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T, F);
    V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(unname(expvals), "black", T, F);
  }else{
    expvals <- rep(0,length(V(bg)$color)); 
    V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- "#00FFFF";
    V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- "#668B8B"
  }
  
  V(bg)$size <- my.rescale(log(cnt2, base=10), 8, 24); 
  
  # layout
  pos.xy <- layout.auto(bg);
  
  # now create the json object
  bnodes <- vector(mode="list");
  node.sizes <- V(bg)$size;
  node.cols <- V(bg)$color;
  node.colsw <- V(bg)$colorw;
  
  shapes <- rep("circle", length(node.nms));
  hit.inx <- node.nms %in% b.mat[,"source"];
  shapes[hit.inx] <- "gene";
  node.lbls = doEntrez2SymbolMapping(node.nms)
  
  for(i in 1:length(node.sizes)){
    bnodes[[i]] <- list(
      id = node.nms[i],
      label=node.lbls[i], 
      size=node.sizes[i], 
      colorb=node.cols[i],
      colorw=node.colsw[i],
      true_size=node.sizes[i], 
      type=shapes[i],
      exp= unname(expvals[node.nms[i]]),
      posx = pos.xy[i,1],
      posy = pos.xy[i,2]
    );
  }
  
  ppi.comps <- vector(mode="list");
  current.net.nm <<- "enrNet"
  ppi.comps[["enrNet"]] <- bg;
  ppi.comps <<- ppi.comps
  
  bedge.mat <- get.edgelist(bg);
  bedge.mat <- cbind(id=1:nrow(bedge.mat), source=bedge.mat[,1], target=bedge.mat[,2]);
  require(RJSONIO);
  initsbls = doEntrez2SymbolMapping(list.genes)
  names(initsbls) = list.genes
  netData <- list(nodes=nodes, edges=edge.mat, bnodes=bnodes, bedges=bedge.mat, enr=enr.mat, id=rownames(enr.mat), sizes=listSizes, hits=hits.query, genelist=initsbls);
  netName = paste0(netNm, ".json");
  sink(netName);
  cat(toJSON(netData));
  sink();
}

GetListEnrGeneNumber <- function(){
  all.enIDs <- NULL;
  listSizes <- list();
  if(anal.type == "genelist"){
    if(numOfLists > 1){
      newDat <- list();
      tot.count <- 0;
      all.nms <- listNms;
      for(i in 1:length(all.nms)){
        dataNm <- all.nms[i];
        dataSet <- readRDS(dataNm);
        gene.mat <- dataSet$prot.mat;
        
        # convert to entrez
        expr.val <- gene.mat[,1];
        en.ids <- rownames(gene.mat);
        
        names(expr.val) <- en.ids;
        newDat[[dataNm]] <- expr.val;
        names(en.ids) <- doEntrez2SymbolMapping(en.ids)
        all.enIDs <- c(all.enIDs, en.ids);
        listSizes[[i]] <- list(
          name = dataNm,
          label = dataNm,
          size = length(en.ids)
          #val = de.prct[i]
        )
      }
      
    }else{
      
      all.enIDs <- rownames(dataSet$prot.mat);
      names(all.enIDs ) <- doEntrez2SymbolMapping(all.enIDs)
      listSizes[[1]] = list(
        name = "datalist1",
        label = "datalist1",
        size = length(all.enIDs)
        #val = de.prct[i]
      )
    }
  }else if(anal.type == "onedata"){
    all.enIDs <- rownames(dataSet$sig.mat);
    names(all.enIDs) <- doEntrez2SymbolMapping(all.enIDs)
    listSizes[[1]] = list(
      name = "dataSet1",
      label = "dataSet1",
      size = length(all.enIDs)
      #val = de.prct[i]
    )
  }else{
    newDat <- list();
    tot.count <- 0;
    listSizes <- list();
    all.nms <- names(mdata.all);
    for(i in 1:length(all.nms)){
      dataNm <- all.nms[i];
      dataSet <- readRDS(dataNm);
      gene.mat <- dataSet$sig.mat;
      
      # convert to entrez
      expr.val <- gene.mat[,1];
      en.ids <- rownames(gene.mat);
      
      names(expr.val) <- en.ids;
      newDat[[dataNm]] <- expr.val;
      names(en.ids) <- doEntrez2SymbolMapping(en.ids)
      all.enIDs <- c(all.enIDs, en.ids);
      listSizes[[i]] <- list(
        name = dataNm,
        label = dataNm,
        size = length(en.ids)
      )
    }
  }
  list.genes <<- all.enIDs;
  listSizes <<- listSizes;
}

InitListEnrichment <- function(type){
  GetListEnrGeneNumber();
  res <- PerformEnrichAnalysis(paste0("enrichment_", type), type, list.genes);
  PrepareEnrichNet(paste0('enrichNet_', type), 'list', "mixed");
  return(res)
}

PerformListEnrichmentView <- function(file.nm, fun.type, netNm, IDs){
  gene.vec <- unlist(strsplit(IDs, "; "));
  gene.vec <- unique(gene.vec);
  sym.vec <- doEntrez2SymbolMapping(gene.vec);
  names(gene.vec) <- sym.vec;
  list.genes <<- gene.vec
  res <- PerformEnrichAnalysis(file.nm, fun.type, list.genes);
  PrepareEnrichNet(netNm, 'list', "mixed");
  return(res);
}

overlap_ratio <- function(x, y, type) {
  x <- unlist(x)
  y <- unlist(y)
  if(type == "mixed"){
    res = 0.5 * length(intersect(x, y))/length(unique(y)) + 0.5 * length(intersect(x, y))/length(unique(c(x,y)))
  }else if(type == "overlap"){
    if(length(x)>length(y)){
      res=length(intersect(x, y))/length(unique(y))
    }else{
      res=length(intersect(x, y))/length(unique(x))
    }
  }else{
    res=length(intersect(x, y))/length(unique(c(x,y)))
  }
  return(res)
}

color_scale <- function(c1="grey", c2="red") {
  pal <- colorRampPalette(c(c1, c2))
  colors <- pal(100)
  return(colors)
}


CalculateDEgeneSetEnr <- function(nms, operation, refNm, filenm){
  nms <- strsplit(nms, ";")[[1]];
  if(anal.type == "metadata" || anal.type == "onedata"){
    com.smbls <- PerformSetOperation_DataEnr(nms, operation, refNm);
  }else{
    com.smbls <- PerformSetOperation_ListEnr(nms, operation, refNm);
  }
  
  sink(filenm);
  cat(toJSON(com.smbls));
  sink();
}

PerformSetOperation_ListEnr <- function(nms, operation, refNm){
  all.nms <- names(mdata.all);
  include.inx <- all.nms %in% nms;
  my.vec <- all.nms[include.inx];
  if(anal.type == "onedata"){
    my.vec = c("1");
  }
  if(operation == "diff"){ # make sure reference is the first
    inx <- which(my.vec == refNm);
    my.vec <- my.vec[-inx];
  }
  com.ids <- NULL;
  ids.list <- list()
  for(i in 1:length(my.vec)){
    if(anal.type != "onedata"){
      dataSet <- readRDS(my.vec[i]);
    }
    if(operation == "diff"){
      ids.list[[i]]=dataSet$GeneAnotDB[,"gene_id"]
      #com.ids <- setdiff(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
    }else if(is.null(com.ids)){
      com.ids <- dataSet$GeneAnotDB[,"gene_id"];
    }else{
      if(operation == "intersect"){
        com.ids <- intersect(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
      }else if(operation == "union"){
        com.ids <- union(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
      }
    }
  }
  if(operation == "diff"){
    dataSet <- readRDS(refNm);
    ids <- unique(unlist(ids.list));
    com.ids <-setdiff(dataSet$GeneAnotDB[,"gene_id"], ids);
  }
  
  com.ids <- unique(as.character(com.ids[!is.na(com.ids)])); # make sure it is interpreted as name not index
  com.symbols <- doEntrez2SymbolMapping(com.ids);
  names(com.symbols) = com.ids
  
  com.symbols<-com.symbols[!is.null(com.symbols)];
  venn.genes <<- com.ids;
  return(com.symbols);
}

PerformSetOperation_DataEnr <- function(nms, operation, refNm){
  
  my.vec <- nms
  if(operation == "diff"){ # make sure reference is the first
    inx <- which(my.vec == refNm);
    my.vec <- my.vec[-inx];
  }
  com.ids <- NULL;
  ids.list <- list()
  if(anal.type == "onedata"){
    my.vec = "dat"
  }
  for(nm in my.vec){
    if(anal.type != "onedata"){
      dataSet <- readRDS(nm);
    }
    if(operation == "diff"){
      ids.list[[nm]]=rownames(dataSet$sig.mat);
      #com.ids <- setdiff(com.ids, dataSet$GeneAnotDB[,"gene_id"]);
    }else if(is.null(com.ids)){
      com.ids <- rownames(dataSet$sig.mat);
    }else{
      if(operation == "intersect"){
        com.ids <- intersect(com.ids, rownames(dataSet$sig.mat));
      }else if(operation=="union"){
        com.ids <- union(com.ids, rownames(dataSet$sig.mat));
      }
    }
  }
  if(operation == "diff"){
    dataSet <- readRDS(refNm);
    ids <- unique(unlist(ids.list));
    com.ids <-setdiff(rownames(dataSet$sig.mat), ids);
  } 
  com.symbols <- doEntrez2SymbolMapping(com.ids);
  names(com.symbols) = com.ids;
  venn.genes <<- com.ids;
  return(com.symbols);
}

PrepareSignatureOfNetworkAnalyst <- function(){
  
  if(anal.type == "genelist"){
    signature.gene <- dataSet$sig.mat;
    
    signature.gene.org <- data.org;
    save(signature.gene, signature.gene.org, file="RShare_networkanalyst.RData");  
  }else if(anal.type == "onedata"){
    if(!file.exists("ExpressResT.rda")){
      return("-1");
    }
    resT <- readRDS("ExpressResT.rda");
    if(exists("P.Value", where=resT)){
      signature.gene <- as.matrix(resT$P.Value);
    }else if(exists("PValue", where=resT)){
      signature.gene <- as.matrix(resT$PValue);
    }
    rownames(signature.gene) <- rownames(resT);
    
    signature.gene.org <- data.org;
    save(signature.gene, signature.gene.org, file="RShare_networkanalyst.RData");  
  }
  
  return(1);
}

##################################################
## R scripts for NetworkAnalyst 
## Description: biological network analysis methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# table.nm is the org code used for sqlite table (ppi)
# for chem type, table.nm is drugbank or ctd
# note, last two param only for STRING database
SearchNetDB <- function(db.type, table.nm, require.exp=TRUE, min.score = 900){ 
  db.typeu <<- db.type
  result.list <- .prepareSigProteinJSON();
  protein.vec <- result.list$protein.vec; # this actually is entrez IDs?
  seed.proteins <<- protein.vec;
  require(RJSONIO);
  # now do the database search
  if(db.type == "ppi"){
    seed.table <- doPpiIDMapping(protein.vec);
    res <- QueryPpiSQLite(table.nm, seed.table$accession, require.exp, min.score);
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    protein.vec <- seed.table$accession;
    edge.res <- data.frame(Source=res[,1],Target=res[,2]);
    row.names(edge.res) <- res[,5];
    write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    node.ids <- c(res[,1], res[,2])
    node.nms <- c(res[,3], res[,4]);
  }else if(db.type == "tf"){ 
    if(table.nm == "encode"){
      table.nm <- paste(table.nm, data.org, sep="_");
    }else{
      table.nm <- toupper(table.nm);
    }
    res <- QueryTFSQLite(table.nm, protein.vec);
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(Source=res[,"entrez"],Target=res[,"tfid"]);
    row.names(edge.res) <- 1:nrow(res);
    write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    node.ids <- c(res[,"entrez"], res[,"tfid"])
    node.nms <- c(res[,"symbol"], res[,"tfname"]);
  }else if(db.type == "mir"){ # in miRNA, table name is org code, colname is id type
    res <- QueryMirSQLite(data.org, "entrez", protein.vec);
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(Source=res[,"entrez"],Target=res[,"mir_acc"]);
    row.names(edge.res) <- 1:nrow(res);
    write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    node.ids <- c(res[,"entrez"], res[,"mir_acc"])
    node.nms <- c(res[,"symbol"], res[,"mir_id"]);
  }else if(db.type == "drug"){
    # note, all drug data is on human, 
    protein.vec <- doEntrez2UniprotMapping(protein.vec);
    protein.vec <- protein.vec[!is.na(protein.vec)];
    res <- QueryDrugSQLite(protein.vec);
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(Source=res[,"upid"],Target=res[,"dbid"]);
    row.names(edge.res) <- 1:nrow(res);
    write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    
    node.ids <- c(res[,"upid"], res[,"dbid"])
    node.nms <- c(res[,"symbol"], res[,"dbname"]);
  }else if(db.type == "disease"){
    # note, all drug data is on human, 
    res <- QueryDiseaseSQLite(protein.vec);
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(Source=res[,"entrez"],Target=res[,"diseaseId"]);
    row.names(edge.res) <- 1:nrow(res);
    write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    
    node.ids <- c(res[,"entrez"], res[,"diseaseId"])
    node.nms <- c(res[,"symbol"], res[,"diseaseName"]);
  }else if(db.type == "tfmir"){
    # note, all drug data is on human, 
    res <- QueryTfmirSQLite(protein.vec);
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(Source=res[,"id1"],Target=res[,"id2"]);
    row.names(edge.res) <- 1:nrow(res);
    write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    
    node.ids <- c(res[,"id1"], res[,"id2"])
    node.nms <- c(res[,"name1"], res[,"name2"]);
  }else if(db.type == "cellcoex"){
    # note, all drug data is on human, 
    res <- QueryCellCoexSQLite(protein.vec);
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(Source=res[,"id1"],Target=res[,"id2"]);
    row.names(edge.res) <- 1:nrow(res);
    write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    
    node.ids <- c(res[,"id1"], res[,"id2"])
    node.nms <- c(res[,"name1"], res[,"name2"]);
  }else if(db.type == "tissueppi"){
    # note, all drug data is on human, 
    res <- QueryDiffNetSQLite(protein.vec);
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(Source=res[,"id1"],Target=res[,"id2"]);
    row.names(edge.res) <- 1:nrow(res);
    write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    
    node.ids <- c(res[,"id1"], res[,"id2"])
    node.nms <- c(res[,"name1"], res[,"name2"]);
  }else if(db.type == "tissuecoex"){
    # note, all drug data is on human, 
    res <- QueryTissueCoexSQLite(protein.vec);
    # no hits
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(Source=res[,"id1"],Target=res[,"id2"]);
    row.names(edge.res) <- 1:nrow(res);
    write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    
    node.ids <- c(res[,"id1"], res[,"id2"])
    node.nms <- c(res[,"name1"], res[,"name2"]);
  }else if(db.type == "chem"){
    res <- QueryChemSQLite(data.org, protein.vec);
    if(nrow(res)==0){ return(c(0,0)); }
    edge.res <- data.frame(Source=res[,"entrez"],Target=res[,"ctdid"]);
    row.names(edge.res) <- 1:nrow(res);
    write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    
    node.ids <- c(res[,"entrez"], res[,"ctdid"])
    node.nms <- c(res[,"symbol"], res[,"name"]);
  }
  
  node.res <- data.frame(Id=node.ids, Label=node.nms);
  node.res <- node.res[!duplicated(node.res$Id),];
  nodeListu <<- node.res
  write.csv(node.res, file="orig_node_list.csv", row.names=FALSE);
  
  ppi.net <<- list(
    db.type=db.type,
    order=1, 
    seeds=protein.vec, 
    table.nm=table.nm, 
    node.data = node.res, 
    edge.data = edge.res,
    require.exp = require.exp,
    min.score = min.score
  );
  
  return(c(nrow(node.res), nrow(res)));
}

.prepareSigProteinJSON <- function(){
  if(anal.type == "genelist"){
    result.list <- .prepareListSeeds();
  }else{ # single expression data or meta.mat
    result.list <- .prepareExpressSeeds();
  }
  return(result.list);
}

.prepareListSeeds <- function(){
  
  protein.list <- list();
  gene.list <- list();
  
  if(numOfLists > 1){
    if(selectedNetDataset %in% c("intersect","union")){
      dataSet = list();
      dataSet$name = selectedNetDataset
      my.vec <- names(mdata.all);
      com.ids <- NULL;
      list.vec <- list()
      for(i in 1:length(my.vec)){
        datSet <- readRDS(my.vec[i]);
        if(is.null(com.ids)){
          com.ids <- datSet$GeneAnotDB[,"gene_id"];
          prot.mat <- datSet$prot.mat
          list.vec[[i]] = com.ids
        }else{
          if(selectedNetDataset == "intersect"){
            com.ids <- datSet$GeneAnotDB[,"gene_id"];
            list.vec[[i]] = com.ids
            #com.ids <- intersect(com.ids, datSet$GeneAnotDB[,"gene_id"]);
          }else{
            com.ids <- union(com.ids, datSet$GeneAnotDB[,"gene_id"]);
          }
          prot.mat <- rbind(prot.mat, as.matrix(datSet$prot.mat[rownames(datSet$prot.mat) %in% com.ids,]))
        }
      }
      if(selectedNetDataset == "intersect"){
        com.ids = Reduce(intersect, list.vec)
        prot.mat <- as.matrix(datSet$prot.mat[rownames(datSet$prot.mat) %in% com.ids,])
      }else{
        com.ids <- unique(as.character(com.ids[!is.na(com.ids)])); # make sure it is interpreted as name not index
      }
      
      com.symbols <- doEntrez2SymbolMapping(com.ids);
      dataSet$GeneAnotDB = data.frame(gene_id=com.ids, accession=com.symbols);
      dataSet$prot.mat = prot.mat;
      dataSet$sig.mat = prot.mat
      dataSet$seeds.proteins = com.ids
    }else{
      my.vec <- names(mdata.all); 
      # make sure reference is the first
      inx <- which(my.vec == selectedNetDataset);
      my.vec <- my.vec[-inx];
      com.ids <- NULL;
      ids.list <- list()
      for(i in 1:length(my.vec)){
        dataSet <- readRDS(my.vec[i]);
        ids.list[[i]]=dataSet$GeneAnotDB[,"gene_id"]
      }
      dataSet <- readRDS(selectedNetDataset);
      ids <- unique(unlist(ids.list));
      com.ids <-setdiff(dataSet$GeneAnotDB[,"gene_id"], ids);
      prot.mat <- as.matrix(dataSet$prot.mat[which(rownames(dataSet$prot.mat) %in% com.ids),])
      com.symbols <- doEntrez2SymbolMapping(com.ids);
      dataSet$GeneAnotDB = data.frame(gene_id=com.ids, accession=com.symbols);
      dataSet$prot.mat = prot.mat;
      dataSet$sig.mat = prot.mat
      dataSet$seeds.proteins = com.ids
    }
  }
  
  # return a json array object
  # each object for a single dataset its sig proteins
  meta.vec <- meta.gene.vec <- meta.seed.expr <- NULL;
  file.create("seed_proteins.txt");
  GeneAnotDB <- NULL;
  
  gene.mat <- dataSet$sig.mat;
  prot.mat <- dataSet$prot.mat;
  write(paste("#DataSet:", dataSet$name),file="sig_genes.txt",append=TRUE);
  write.table(dataSet$sig.mat, file="sig_genes.txt", append=TRUE);
  
  meta.gene.vec <- c(meta.gene.vec, rownames(gene.mat));
  gene.list[[dataSet$name]] <- list(gene=rownames(gene.mat),logFC=unname(gene.mat[,1]));
  GeneAnotDB <- rbind(GeneAnotDB, dataSet$GeneAnotDB);
  meta.seed.expr <- c(meta.seed.expr, prot.mat[,1]);
  write(paste("#DataSet:", dataSet$name),file="seed_proteins.txt",append=TRUE);
  write.table(cbind(Emblprotein=rownames(prot.mat), Expression=prot.mat[,1]), file="seed_proteins.txt", row.names=F, quote=F,append=TRUE);
  protein.vec <- prot.mat[,1];
  meta.vec <- c(meta.vec, names(protein.vec));
  if(length(protein.vec) == 1){
    protein.vec <- as.matrix(protein.vec)
  }   
  protein.list[[dataSet$name]] <- signif(protein.vec, 3);
  
  gene.list$name <- dataSet$name;
  seed.genes <<- unique(meta.gene.vec);
  
  meta.seed.df <- as.matrix(meta.seed.expr);
  rownames(meta.seed.df) <- names(meta.seed.expr);
  
  seed.expr <- RemoveDuplicates(meta.seed.df, "max", quiet=F); 
  seed.expr <<- seed.expr[,1];
  protein.vec <- unique(meta.vec);
  
  result = list(
    gene.list = gene.list,
    protein.list = protein.list,
    protein.vec = protein.vec
  );
  return(result)
}


# 2nd order network, only for PPI
ExpandNetworkSearch <- function(){
  if(ppi.net$order == 0){
    SearchNetDB(ppi.net$db.type, ppi.net$table.nm);
    CreateGraph();
  }
  protein.vec <- V(overall.graph)$name;
  table.nm <- ppi.net$table.nm;
  if(db.typeu == "ppi"){
    res <- QueryPpiSQLite(table.nm, protein.vec, ppi.net$require.exp, ppi.net$min.score);
  }else if (db.typeu == "tissuecoex"){
    res <- QueryTissueCoexSQLite(protein.vec);
  }else if (db.typeu == "tissueppi"){
    res <- QueryDiffNetSQLite(protein.vec);
  }else if (db.typeu == "cellcoex"){
    res <- QueryCellCoexSQLite(protein.vec);
  }
  edge.res <- data.frame(Source=res[,1],Target=res[,2]);
  #row.names(edge.res) <- res[,5];
  
  node.ids <- c(res[,1], res[,2])
  node.nms <- c(res[,3], res[,4]);
  node.res <- data.frame(Id=node.ids, Label=node.nms);
  node.res <- node.res[!duplicated(node.res$Id),];
  if(nrow(node.res) < 10000){ # only overwrite if within the range
    write.csv(edge.res, file="orig_edge_list.csv",row.names=FALSE);
    write.csv(node.res, file="orig_node_list.csv", row.names=FALSE);
    ppi.net <<- list(db.type=db.typeu, order=2, seeds=protein.vec, table.nm=table.nm, node.data = node.res, edge.data = edge.res);
  }
  return(nrow(node.res));
}

# zero-order network - create ppi nets from only input (seeds)
BuildSeedProteinNet <- function(){
  
  nodes = V(overall.graph)$name;
  hit.inx <-  nodes %in% seed.proteins;
  nodes2rm <- nodes[!hit.inx];
  g <- simplify(delete.vertices(overall.graph, nodes2rm));
  
  nodeList <- get.data.frame(g, "vertices");
  nodeList <- nodeList[,1:2];
  colnames(nodeList) <- c("Id", "Label");
  
  write.csv(nodeList, file="orig_node_list.csv", quote=F);
  nd.inx <- ppi.net$node.data[,1] %in% nodeList[,1];
  
  edgeList <- get.data.frame(g, "edges");
  edgeList <- edgeList[,1:2];
  colnames(edgeList) <- c("Source", "Target");
  write.csv(edgeList, file="orig_edge_list.csv",  quote=F);
  
  # update ppi.net 
  ppi.net$order = 0;
  ppi.net$node.data <- nodeList;
  ppi.net$edge.data <- edgeList;
  ppi.net <<- ppi.net;
}

# create igraph from the edgelist saved from graph DB
# and decompose into subnets
CreateGraph <- function(){
  
  require('igraph');
  node.list <- ppi.net$node.data;
  edge.list <- ppi.net$edge.data;
  
  seed.proteins <- ppi.net$seeds;
  overall.graph <- simplify(graph.data.frame(edge.list, directed=FALSE, vertices=node.list));
  
  # add node expression value
  if(ppi.net$db.type == "ppi"){
    ## Temp fix seed.expr all in entrez? 
    ## all converted to new IDs based on PPI from the given organism 
    newIDs <- doPpiIDMapping(names(seed.expr))$accession;    
  }else{# all entrez in mirNet
    newIDs <- names(seed.expr);
  }
  match.index <- match(V(overall.graph)$name, newIDs);
  expr.vals <- seed.expr[match.index];
  #   expr.vec <- abs(expr.vals) # logFC can be negative!!
  expr.vec <- expr.vals;
  names(expr.vec)= V(overall.graph)$name;
  expr.vec <<- expr.vec[!is.na(expr.vec)]
  overall.graph <- set.vertex.attribute(overall.graph, "abundance", index = V(overall.graph), value = expr.vals);
  
  hit.inx <- seed.proteins %in% node.list[,1];
  seed.proteins <<- seed.proteins[hit.inx];
  
  substats <- DecomposeGraph(overall.graph);
  overall.graph <<- overall.graph;
  if(!is.null(substats)){
    return(c(length(seed.genes), length(seed.proteins), nrow(node.list), nrow(edge.list), length(ppi.comps), substats));        
  }else{
    return(0);
  }
}

FilterBipartiNet <- function(nd.type, min.dgr, min.btw){
  
  all.nms <- V(overall.graph)$name;
  edge.mat <- get.edgelist(overall.graph);
  dgrs <- degree(overall.graph);
  nodes2rm.dgr <- nodes2rm.btw <- NULL;
  
  if(nd.type == "gene"){
    hit.inx <- all.nms %in% edge.mat[,1];
  }else if(nd.type=="other"){
    hit.inx <- all.nms %in% edge.mat[,2];
  }else{ # all
    hit.inx <- rep(TRUE, length(all.nms));
  }
  
  if(min.dgr > 0){
    rm.inx <- dgrs <= min.dgr & hit.inx;
    nodes2rm.dgr <- V(overall.graph)$name[rm.inx];
  }
  if(min.btw > 0){
    btws <- betweenness(overall.graph);
    rm.inx <- btws <= min.btw & hit.inx;
    nodes2rm.btw <- V(overall.graph)$name[rm.inx];
  }
  
  nodes2rm <- unique(c(nodes2rm.dgr, nodes2rm.btw));
  overall.graph <- simplify(delete.vertices(overall.graph, nodes2rm));
  current.msg <<- paste("A total of", length(nodes2rm) , "was reduced.");
  substats <- DecomposeGraph(overall.graph);
  if(!is.null(substats)){
    overall.graph <<- overall.graph;
    return(c(length(seed.genes),length(seed.proteins), vcount(overall.graph), ecount(overall.graph), length(ppi.comps), substats));
  }else{
    return(0);
  }
}

PrepareNetwork <- function(net.nm, json.nm){
  
  my.ppi <- ppi.comps[[net.nm]];
  nd.nms <- V(my.ppi)$name;
  
  
  convertIgraph2JSON(net.nm, json.nm);
  current.net.nm <<- net.nm;
  return(1);
}

PrepareNetworkUpload <- function(net.nm, json.nm){
  
  my.ppi <- ppi.comps[[net.nm]];
  nd.nms <- V(my.ppi)$name;
  
  convertIgraph2JSONFromFile(net.nm, json.nm);
  current.net.nm <<- net.nm;
  return(1);
}

GetNodeIDs <- function(){
  V(overall.graph)$name;
}

GetNodeNames <- function(){
  V(overall.graph)$Label;
}

GetNodeDegrees <- function(){
  degree(overall.graph);
}

GetNodeBetweenness <- function(){
  round(betweenness(overall.graph, directed=F, normalized=F), 2);
}

# use microservice
GetPCSFNet <- function(){
  
  print("Peforming PCSF ....");
  
  library(RSclient);
  rsc <- RS.connect();
  RS.assign(rsc, "my.dir", getwd()); 
  RS.eval(rsc, setwd(my.dir));
  
  dat.out <- list(data=overall.graph, terminals = expr.vec);
  RS.assign(rsc, "dat.in", dat.out); 
  my.fun <- function(){
    require('PCSF');
    edg <- get.edgelist(dat.in$data);
    edg <- as.data.frame(edg);
    edg$V3 <- rep(1, nrow(edg));
    colnames(edg) <- c("from", "to", "cost");
    ppi <- construct_interactome(edg);
    g <- PCSF(ppi, dat.in$terminals, w = 5, b = 100, mu = 0.0005);
    return(g);
  }
  
  RS.assign(rsc, my.fun);
  g <-  RS.eval(rsc, my.fun());
  RS.close(rsc);
  
  nodeList <- get.data.frame(g, "vertices");
  colnames(nodeList) <- c("Id", "Label");
  write.csv(nodeList, file="orig_node_list.csv", row.names=F, quote=F);
  
  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  write.csv(edgeList, file="orig_edge_list.csv", row.names=F, quote=F);
  
  path.list <- NULL;
  substats <- DecomposeGraph(g);
  if(!is.null(substats)){
    overall.graph <<- g;
    return(c(length(seed.genes),length(seed.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), substats));
  }else{
    return(0);
  }
}

# for a given graph, obtain the smallest subgraphs that contain
# all the seed nodes. This is acheived by iteratively remove 
# the marginal nodes (degree = 1) that are not in the seeds
GetMinConnectedGraphs <- function(max.len = 200){
  set.seed(8574);
  # first get shortest paths for all pair-wise seeds
  my.seeds <- seed.proteins;
  sd.len <- length(my.seeds);
  paths.list <-list();
  
  # first trim overall.graph to remove no-seed nodes of degree 1
  dgrs <- degree(overall.graph);
  keep.inx <- dgrs > 1 | (names(dgrs) %in% my.seeds);
  nodes2rm <- V(overall.graph)$name[!keep.inx];
  overall.graph <-  simplify(delete.vertices(overall.graph, nodes2rm));
  
  # need to restrict the operation b/c get.shortest.paths is very time consuming
  # for top max.len highest degrees
  if(sd.len > max.len){
    hit.inx <- names(dgrs) %in% my.seeds;
    sd.dgrs <- dgrs[hit.inx];
    sd.dgrs <- rev(sort(sd.dgrs));
    # need to synchronize all (seed.proteins) and top seeds (my.seeds)
    seed.proteins <- names(sd.dgrs);    
    my.seeds <- seed.proteins[1:max.len];
    sd.len <- max.len;
    current.msg <<- paste("The minimum connected network was computed using the top", sd.len, "seed proteins in the network based on their degrees.");
  }else{
    current.msg <<- paste("The minimum connected network was computed using all seed proteins in the network.");
  }
  # now calculate the shortest paths for 
  # each seed vs. all other seeds (note, to remove pairs already calculated previously)
  for(pos in 1:sd.len){
    paths.list[[pos]] <- get.shortest.paths(overall.graph, my.seeds[pos], seed.proteins[-(1:pos)])$vpath;
  }
  nds.inxs <- unique(unlist(paths.list));
  nodes2rm <- V(overall.graph)$name[-nds.inxs];
  g <- simplify(delete.vertices(overall.graph, nodes2rm));
  
  nodeList <- get.data.frame(g, "vertices");
  colnames(nodeList) <- c("Id", "Label");
  write.csv(nodeList, file="orig_node_list.csv", row.names=F, quote=F);
  
  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  write.csv(edgeList, file="orig_edge_list.csv", row.names=F, quote=F);
  
  path.list <- NULL;
  substats <- DecomposeGraph(g);
  if(!is.null(substats)){
    overall.graph <<- g;
    return(c(length(seed.genes),length(seed.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), substats));
  }else{
    return(0);
  }
}

UpdateSubnetStats <- function(){
  old.nms <- names(ppi.comps);
  net.stats <- ComputeSubnetStats(ppi.comps);
  ord.inx <- order(net.stats[,1], decreasing=TRUE);
  net.stats <- net.stats[ord.inx,];
  rownames(net.stats) <- old.nms[ord.inx];
  net.stats <<- net.stats;
}

# exclude nodes in current.net (networkview)
ExcludeNodes <- function(nodeids, filenm){
  
  nodes2rm <- strsplit(nodeids, ";")[[1]];
  current.net <- ppi.comps[[current.net.nm]];
  current.net <- delete.vertices(current.net, nodes2rm);
  
  # need to remove all orphan nodes
  bad.vs<-V(current.net)$name[degree(current.net) == 0];
  current.net <- delete.vertices(current.net, bad.vs);
  
  # return all those nodes that are removed 
  nds2rm <- paste(c(bad.vs, nodes2rm), collapse="||");
  
  # update topo measures
  node.btw <- as.numeric(betweenness(current.net));
  node.dgr <- as.numeric(degree(current.net));
  node.exp <- as.numeric(get.vertex.attribute(current.net, name="abundance", index = V(current.net)));
  nms <- V(current.net)$name;
  hit.inx <- match(nms, ppi.net$node.data[,1]);
  lbls <- ppi.net$node.data[hit.inx,2];
  
  nodes <- vector(mode="list");
  for(i in 1:length(nms)){
    nodes[[i]] <- list(
      id=nms[i], 
      label=lbls[i],
      degree=node.dgr[i], 
      between=node.btw[i],
      expr = node.exp[i]
    );
  }
  # now only save the node pos to json
  require(RJSONIO);
  netData <- list(deletes=nds2rm,nodes=nodes);
  sink(filenm);
  cat(toJSON(netData));
  sink();
  
  ppi.comps[[current.net.nm]] <<- current.net;
  UpdateSubnetStats();
  
  # remember to forget the cached layout, and restart caching, as this is now different object (with the same name)
  #forget(PerformLayOut_mem);
  return(filenm);
}

# exclude nodes in overall net (network builder)
ExcludeNodesOverall <- function(nodeids, id.type, vismode){
  # all convert to uniprot ID 
  lines <- strsplit(nodeids, "\r|\n|\r\n")[[1]];
  lines<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", lines, perl=TRUE);
  if(vismode != "network"){
    prot.anots <- convertIdToEntrez(lines, id.type);
    nodes2rm <- unique(prot.anots$accession);
  }else{
    prot.anots <- lines
    nodes2rm <- unique(lines);
  }
  
  # now find the overlap
  nodes2rm <- nodes2rm[nodes2rm %in% V(overall.graph)$name];
  g <- delete.vertices(overall.graph, nodes2rm);
  
  nodeList <- get.data.frame(g, "vertices");
  colnames(nodeList) <- c("Id", "Label");
  
  write.csv(nodeList, file="orig_node_list.csv", row.names=F, quote=F);
  
  edgeList <- get.data.frame(g, "edges");
  edgeList <- cbind(rownames(edgeList), edgeList);
  colnames(edgeList) <- c("Id", "Source", "Target");
  write.csv(edgeList, file="orig_edge_list.csv", row.names=F, quote=F);
  
  substats <- DecomposeGraph(g);
  if(!is.null(substats)){
    overall.graph <<- g;
    #forget(PerformLayOut_mem);
    return(c(length(seed.genes),length(seed.proteins), nrow(nodeList), nrow(edgeList), length(ppi.comps), substats));
  }else{
    return(0);
  }
}

PrepareSubnetDownloads <- function(nm){
  g <- ppi.comps[[nm]];
  # need to update graph so that id is compound names rather than ID
  V(g)$name <- as.character(doID2LabelMapping(V(g)$name));
  saveNetworkInSIF(g, nm);
}

# adapted from BioNet
saveNetworkInSIF <- function(network, name){
  edges <- .graph.sif(network=network, file=name);
  sif.nm <- paste(name, ".sif", sep="");
  if(length(list.edge.attributes(network))!=0){
    edge.nms <- .graph.eda(network=network, file=name, edgelist.names=edges);
    sif.nm <- c(sif.nm, edge.nms);
    
  }
  if(length(list.vertex.attributes(network))!=0){
    node.nms <- .graph.noa(network=network, file=name);
    sif.nm <- c(sif.nm, node.nms);
  }
  # need to save all sif and associated attribute files into a zip file for download
  zip(paste(name,"_sif",".zip", sep=""), sif.nm);
}

# internal function to write cytoscape .sif file
.graph.sif <- function(network, file){
  edgelist.names <- igraph::get.edgelist(network, names=TRUE)
  edgelist.names <- cbind(edgelist.names[,1], rep("pp", length(E(network))), edgelist.names[,2]);
  write.table(edgelist.names, row.names=FALSE, col.names=FALSE, file=paste(file, ".sif", sep=""), sep="\t", quote=FALSE)
  return(edgelist.names) 
}

doID2LabelMapping <- function(entrez.vec){
  if(exists("nodeListu")){
    hit.inx <- match(entrez.vec, nodeListu[, "Id"]);
    symbols <- nodeListu[hit.inx, "Label"];
    
    # if not gene symbol, use id by itself
    na.inx <- is.na(symbols);
    symbols[na.inx] <- entrez.vec[na.inx];
    return(symbols);
  }else{ # network upload
    return(entrez.vec);
  }
} 

# internal method to write cytoscape node attribute files
.graph.noa <- function(network, file){
  all.nms <- c();
  attrib <- list.vertex.attributes(network)
  for(i in 1:length(attrib)){
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "character")
    {
      type <- "String"
    }
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "integer")
    {
      type <- "Integer"
    }
    if(is(get.vertex.attribute(network, attrib[i]))[1] == "numeric")
    {
      type <- "Double"
    }
    noa <- cbind(V(network)$name, rep("=", length(V(network))), get.vertex.attribute(network, attrib[i]))
    first.line <- paste(attrib[i], " (class=java.lang.", type, ")", sep="")
    file.nm <- paste(file, "_", attrib[i], ".NA", sep="");
    write(first.line, file=file.nm, ncolumns = 1, append=FALSE, sep=" ")
    write.table(noa, row.names = FALSE, col.names = FALSE, file=file.nm, sep=" ", append=TRUE, quote=FALSE);
    all.nms <- c(all.nms, file.nm);
  }
  return(all.nms);
}

# internal method to write cytoscape edge attribute files
.graph.eda <- function(network, file, edgelist.names){
  all.nms <- c();
  attrib <- list.edge.attributes(network)
  for(i in 1:length(attrib)){
    if(is(get.edge.attribute(network, attrib[i]))[1] == "character")
    {
      type <- "String"
    }
    if(is(get.edge.attribute(network, attrib[i]))[1] == "integer")
    {
      type <- "Integer"
    }
    if(is(get.edge.attribute(network, attrib[i]))[1] == "numeric")
    {
      type <- "Double"
    }
    eda <- cbind(cbind(edgelist.names[,1], rep("(pp)", length(E(network))), edgelist.names[,3]), rep("=", length(E(network))), get.edge.attribute(network, attrib[i]))
    first.line <- paste(attrib[i], " (class=java.lang.", type, ")", sep="");
    file.nm <- paste(file, "_", attrib[i], ".EA", sep="");
    write(first.line, file=file.nm, ncolumns=1, append=FALSE, sep =" ")
    write.table(eda, row.names = FALSE, col.names = FALSE, file=file.nm, sep=" ", append=TRUE, quote=FALSE);
    all.nms <- c(all.nms, file.nm);
  }
  return(all.nms);
}

SetCellCoexNumber <-function(num){
  cellCoexNumber <<- num;
}

SetDiffNetName <-function(nm){
  diffNetName <<- nm;
}

SetDiffFilter <-function(pct){
  diffPct <<- pct/10;
}

# note: hit.query, resTable must synchronize
PerformNetEnrichment <- function(file.nm, fun.type, IDs){
  # prepare query
  ora.vec <- NULL;
  if(ppi.net$db.type == 'ppi'){
    if(data.org == "ath"){
      idtype <- "tair"
    }else if(data.org %in% c("bsu", "tbr", "cel", "dme", "eco", "pfa", "pae") & net.type == "string"){
      idtype <- "string"
    }else if(data.org %in% c("bta","dre","rno","gga","hsa","mmu") & net.type == "string"){
      idtype <- "emblprotein"
    }else if(data.org %in% c("hsa","mmu", "cel", "dme","sce") & net.type %in% c("innate", "irefinx", "rolland")){
      idtype <- "uniprot"
    }else if(data.org == "sce" & net.type == "string"){ # only for yeast
      idtype <- "emblgene";
    }
    if(idtype=="uniprot"){
      uniprot.vec <- unlist(strsplit(IDs, "; "));
      ora.vec <- doUniprot2EntrezMapping(uniprot.vec);
      names(ora.vec) <- uniprot.vec;
    }else if(idtype=="emblprotein"){
      emblprotein.vec <- unlist(strsplit(IDs, "; "))
      ora.vec <- doEmblProtein2EntrezMapping(emblprotein.vec);
      names(ora.vec) <- emblprotein.vec;
    }else if(idtype=="string"){
      string.vec <- unlist(strsplit(IDs, "; "))
      ora.vec <- doString2EntrezMapping(string.vec);
      names(ora.vec) <- string.vec;
    }else if(idtype=="emblgene"){
      emblgene.vec <- unlist(strsplit(IDs, "; "))
      ora.vec <- doEmblGene2EntrezMapping(emblgene.vec);
      names(ora.vec) <- emblgene.vec;
    }else{
      ora.vec <- unlist(strsplit(IDs, "; "));
      names(ora.vec) <- ora.vec;
    }
  }else{ # net is tf/mir/drug, they already in entrez
    ora.vec <- unlist(strsplit(IDs, "; "));
    names(ora.vec) <- as.character(ora.vec);
  }
  res <- PerformEnrichAnalysis(file.nm, fun.type, ora.vec);
  return(res);
  
}

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
  tf.db <- dbConnect(SQLite(), paste(sqlite.path, "tfmir.sqlite", sep="")); 
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
  chem.db <- dbConnect(SQLite(), paste(sqlite.path, "tfac.sqlite", sep="")); 
  query <- paste (shQuote(q.vec),collapse=",");
  statement <- paste("SELECT * FROM ", table.nm, " WHERE entrez IN (",query,")", sep="");
  return(.query.sqlite(chem.db, statement));
}

doPpiIDMapping <- function(q.vec){
  if(data.org == "ath"){
    db.map <-  queryGeneDB("tair", data.org);
  }else if(data.org == "sce"){
    if(net.type == "string"){ # only for yeast
      db.map <-  queryGeneDB("entrez_embl_gene", data.org);
    }else{
      db.map <-  queryGeneDB("entrez_uniprot", data.org);
    }
  }else if(net.type %in% c("innate", "irefinx", "rolland")){
    db.map <-  queryGeneDB("entrez_uniprot", data.org);
  }else{
    if(data.org %in% c("bsu", "tbr", "cel", "dme", "eco", "pfa","pae")){
      db.map <-  queryGeneDB("entrez_string", data.org);
    }else if (data.org %in% c("mmu","hsa") && net.type == "string"){
      db.map <-  queryGeneDB("entrez_string", data.org);
    }else if(data.org == "mtb"){
      db.map <-  queryGeneDB("entrez", data.org);
    }else{
      db.map <-  queryGeneDB("entrez_embl_protein", data.org);
    }
  }
  hit.inx <- match(q.vec, db.map[, "gene_id"]);
  ppi.mat <- db.map[hit.inx, ]; 
  
  # fix the factor col related to library issue
  i <- sapply(ppi.mat, is.factor)
  ppi.mat[i] <- lapply(ppi.mat[i], as.character)
  if(data.org %in% c("pae", "mtb")){
    ppi.mat = ppi.mat[,c(2,1)];
    colnames(ppi.mat) = c("gene_id", "accession");
  }
  return(ppi.mat);
}

doUniprot2EntrezMapping <- function(uniprot.vec){
  # db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep="");
  # db.map <-  readRDS(db.path);
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

##################################################
## R script for NetworkAnalyst
## Description: prepare data for Venn diagram
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# create a list store all possible combination (for a max of 4)
# note, the whole region are divided into small areas (16 for 4; 7 for 3, 3 for 2)
# for instance:
# a: a unique (no b, no c)
# ab: a and b, no c
PrepareVennData<-function(){
  newDat <- list();
  if(anal.type == "metadata"){
    sel.inx <- mdata.all==1;
    sel.nms <- names(mdata.all)[sel.inx];
  }else{
    sel.nms <- listNms;
  }
  venn.genenb <- vector()
  sel.dats <- list();
  for(i in 1:length(sel.nms)){
    nm = sel.nms[i]
    dataSet <- readRDS(nm);
    if(anal.type == "metadata"){
      sel.dats[[nm]] <- rownames(dataSet$sig.mat)
    }else{
      sel.dats[[nm]] <- rownames(dataSet$prot.mat)
    }
    venn.genenb[i] = length(sel.dats[[nm]])
  }
  if(anal.type == "metadata" & meta.selected){
    sel.dats[["meta_dat"]] <- as.character(meta.stat$de);
    venn.genenb[length(venn.genenb) + 1] = length(as.character(meta.stat$de))
  }
  if(length(sel.dats) == 2){
    venn.list <<- Prepare2Venn(sel.dats);
  }else if(length(sel.dats) == 3){
    venn.list <<- Prepare3Venn(sel.dats);
  }else if(length(sel.dats) == 4){
    venn.list <<- Prepare4Venn(sel.dats);
  }else{
    venn.list <<- Prepare4Venn(sel.dats[c(1:4)]);
  }
  venn.list <<- sel.dats;
  venn.genenb <<- venn.genenb
  venn.list.up <<- sel.dats;
  venn.genenb.up <<- venn.genenb
  return(1);
}

PrepareSelVennData<-function(selectedNms){
  newDat <- list();
  sel.nms <- unlist(strsplit(selectedNms, ";"));
  venn.genenb <- vector()
  sel.dats <- list();
  for(i in 1:length(sel.nms)){
    nm = sel.nms[i]
    if(nm != "meta_dat"){
      dataSet <- readRDS(nm);
      if(anal.type == "metadata"){
        sel.dats[[nm]] <- rownames(dataSet$sig.mat)
      }else{
        sel.dats[[nm]] <- rownames(dataSet$prot.mat)
      }
      venn.genenb[i] = length(sel.dats[[nm]])
    }else{
      sel.dats[[nm]] <- as.character(meta.stat$de);
      venn.genenb[i] = length(as.character(meta.stat$de))
    }
  }
  if(length(sel.dats) == 2){
    venn.list <- Prepare2Venn(sel.dats);
  }else if(length(sel.dats) == 3){
    venn.list <- Prepare3Venn(sel.dats);
  }else if(length(sel.dats) == 4){
    venn.list <- Prepare4Venn(sel.dats);
  }
  venn.list.up <<- sel.dats;
  venn.genenb.up <<- venn.genenb
  return(1);
}


#3
Prepare2Venn <- function(dat){
  nms <- names(dat);  
  a <- nms[1];
  b <- nms[2];
  ab <- paste(a, b, sep="");
  
  a.l <- dat[[a]];
  b.l <- dat[[b]];
  
  vennData <- list();
  vennData[[a]] <- setdiff(a.l, b.l);
  vennData[[b]] <- setdiff(b.l, a.l);    
  vennData[[ab]] <- intersect(b.l, a.l);
  vennData <<- vennData;
}

#7
Prepare3Venn <- function(dat){
  nms <- names(dat);  
  a <- nms[1];
  b <- nms[2];
  c <- nms[3];
  ab <- paste(a, b, sep="");
  ac <- paste(a, c, sep="");
  bc <- paste(b, c, sep="");
  abc <- paste(a, b, c, sep="");
  
  a.l <- dat[[a]];
  b.l <- dat[[b]];
  c.l <- dat[[c]];
  
  vennData <- list();
  vennData[[a]] <- setdiff(a.l, union(b.l, c.l));
  vennData[[b]] <- setdiff(b.l, union(a.l, c.l));    
  vennData[[c]] <- setdiff(c.l, union(a.l, b.l));    
  vennData[[ab]] <- setdiff(intersect(a.l, b.l), c.l);
  vennData[[ac]] <- setdiff(intersect(a.l, c.l), b.l);
  vennData[[bc]] <- setdiff(intersect(b.l, c.l), a.l);
  vennData[[abc]] <- intersect(intersect(a.l, b.l), c.l);
  vennData <<- vennData;
}

# 15
Prepare4Venn <- function(dat){
  nms <- names(dat);  
  a <- nms[1];
  b <- nms[2];
  c <- nms[3];
  d <- nms[4];
  ab <- paste(a, b, sep="");
  ac <- paste(a, c, sep="");
  ad <- paste(a, d, sep="");
  bc <- paste(b, c, sep="");
  bd <- paste(b, d, sep="");
  cd <- paste(c, d, sep="");
  abc <- paste(a, b, c, sep="");
  abd <- paste(a, b, d, sep="");
  acd <- paste(a, c, d, sep="");
  bcd <- paste(b, c, d, sep="");
  abcd <- paste(a, b, c, d, sep="");
  
  a.l <- dat[[a]];
  b.l <- dat[[b]];
  c.l <- dat[[c]];
  d.l <- dat[[d]];
  
  vennData <- list();
  vennData[[a]] <- setdiff(a.l, unique(c(b.l, c.l, d.l)));
  vennData[[b]] <- setdiff(b.l, unique(c(a.l, c.l, d.l)));    
  vennData[[c]] <- setdiff(c.l, unique(c(a.l, b.l, d.l)));    
  vennData[[d]] <- setdiff(d.l, unique(c(a.l, b.l, c.l))); 
  vennData[[ab]] <- setdiff(intersect(a.l, b.l), union(c.l, d.l));
  vennData[[ac]] <- setdiff(intersect(a.l, c.l), union(b.l, d.l));
  vennData[[ad]] <- setdiff(intersect(a.l, d.l), union(b.l, c.l));
  vennData[[bc]] <- setdiff(intersect(b.l, c.l), union(a.l, d.l));
  vennData[[bd]] <- setdiff(intersect(b.l, d.l), union(a.l, c.l));
  vennData[[cd]] <- setdiff(intersect(c.l, d.l), union(a.l, b.l));
  vennData[[abc]] <- setdiff(intersect(intersect(a.l, b.l), c.l), d.l);
  vennData[[abd]] <- setdiff(intersect(intersect(a.l, b.l), d.l), c.l);
  vennData[[acd]] <- setdiff(intersect(intersect(a.l, c.l), d.l), b.l);
  vennData[[bcd]] <- setdiff(intersect(intersect(b.l, c.l), d.l), a.l);
  vennData[[abcd]] <- intersect(intersect(a.l, b.l), intersect(c.l, d.l));
  vennData <<- vennData;
}

GetSelectedDataNumber <- function(){
  return(length(venn.list));
}

GetSelectedDataNames <- function(){
  return(paste(names(venn.list), collapse=";"));
}

GetSelectedDataGeneNumber<- function(){
  return(paste(venn.genenb, collapse=";"));
}

GetSelectedDataNumberUpdated <- function(){
  return(length(venn.list.up));
}

GetSelectedDataNamesUpdated <- function(){
  return(paste(names(venn.list.up), collapse=";"));
}

GetSelectedDataGeneNumberUpdated<- function(){
  return(paste(venn.genenb.up, collapse=";"));
}


#areas is allname concated
GetVennGeneNames <- function(areas){
  nms <- strsplit(areas, "\\|\\|")[[1]];
  gene.vec <- NULL;
  for(nm in nms){
    gene.vec <- c(gene.vec, vennData[[nm]]);
  }
  gene.vec <- unique(gene.vec);
  # from entrez to symbols
  sym.vec <- doEntrez2SymbolMapping(gene.vec);
  names(gene.vec) <- sym.vec;
  venn.genes <<- gene.vec;
  return(paste(unique(sym.vec), collapse="||"));
}

PerformVennEnrichment <- function(file.nm, fun.type){
  res <- PerformEnrichAnalysis(file.nm, fun.type, venn.genes);
  return(res);
}