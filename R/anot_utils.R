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
    GeneAnotDB <-doProteinIDMapping(rownames(gene.mat), type);
    
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
  feature.vec <- id.vec;
  if(idType %in% c("entrez", "symbol", "refseq", "genbank", "emblgene","emblprotein", "embltranscript", "orfid", "tair", "wormbase")){
    anot.id <- doGeneIDMapping(feature.vec, idType);
  }else{
    anot.id <- doProbeMapping(feature.vec, idType);
  }
  names(anot.id) <- id.vec;
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

# mapping between genebank, refseq and entrez
doGeneIDMapping <- function(q.vec, type){
  if(is.null(q.vec)){
    db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
    db.map <-  readRDS(db.path);
    q.vec <- db.map[, "gene_id"];
    type = "entrez";
  }
  
  if(type == "symbol"){
    db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
    db.map <-  readRDS(db.path);
    hit.inx <- match(q.vec, db.map[, "symbol"]);
  }else if(type == "entrez"){
    db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
    db.map <-  readRDS(db.path);
    hit.inx <- match(q.vec, db.map[["gene_id"]]);
  }else{
    if(type == "embltranscript"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_transcript.rds", sep="");
    }else {
      # note, some ID can have version number which is not in the database
      # need to strip it off NM_001402.5 => NM_001402
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      if(type == "genbank"){
        db.path <- paste(lib.path, data.org, "/entrez_gb.rds", sep="");
      }else if(type == "refseq"){
        db.path <- paste(lib.path, data.org, "/entrez_refseq.rds", sep="");
      }else if(type == "emblgene"){
        db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
      }else if(type == "embltranscript"){
        db.path <- paste(lib.path, data.org, "/entrez_embl_transcript.rds", sep="");
      }else if(type == "emblprotein"){
        db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
      }else if(type == "orfid"){ # only for yeast
        db.path <- paste(lib.path, data.org, "/entrez_orf.rds", sep="");
      }else if(type == "tair"){ # only for ath
        db.path <- paste(lib.path, data.org, "/tair.rds", sep="");
      }else if(type == "wormbase"){ # only for cel
        db.path <- paste(lib.path, data.org, "/entrez_wormbase.rds", sep="");
      }else{
        print("Unknown data type");
        return(0);
      }
    }
    db.map <-  readRDS(db.path);
    hit.inx <- match(q.vec, db.map[["accession"]]);
  }
  # if(type == "wormbase"){
  #   entrezs=as.vector(db.map[hit.inx, "gene_id"]);
  #   entrezs=entrezs$gene_id
  #   mode(entrezs) <- "character";
  #   rm(db.map, q.vec); gc();
  #   return(entrezs);
  # } else {
  entrezs=db.map[hit.inx, "gene_id"];
  mode(entrezs) <- "character";
  rm(db.map, q.vec); gc();
  return(entrezs);
  # }
}

doEntrez2SymbolMapping <- function(entrez.vec){
  db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
  gene.map <- readRDS(db.path);
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
  db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
  gene.map <- readRDS(db.path);
  
  hit.inx <- match(entrez.vec, gene.map[, "gene_id"]);
  anot.mat <- gene.map[hit.inx, c("gene_id", "symbol", "name")];
  
  na.inx <- is.na(hit.inx);
  anot.mat[na.inx, "symbol"] <- entrez.vec[na.inx];
  anot.mat[na.inx, "name"] <- 'NA';
  return(anot.mat);
}

doProteinIDMapping <- function(q.vec, type){
  if(type == "entrez"){
    # need to get only our data
    db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
    db.map <-  readRDS(db.path);
    hit.inx <- match(q.vec, db.map[, "gene_id"]);
    entrezs <- db.map[hit.inx, ]
    entrezs <- entrezs[c(1,1)]
    colnames(entrezs) = c("gene_id", "accession");
  }else if(type == "symbol"){
    db.path <- paste(lib.path, data.org, "/entrez.rds", sep="");
    gene.map <- readRDS(db.path);
    hit.inx <- match(q.vec, gene.map[, "symbol"]);
    entrezs <- gene.map[hit.inx, ];
    entrezs = entrezs[c(1,2)];
    colnames(entrezs) <- c("gene_id", "accession")     
  }else{
    if(type == "genbank"){
      # note, some ID can have version number which is not in the database
      # need to strip it off NM_001402.5 => NM_001402
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.path <- paste(lib.path, data.org, "/entrez_gb.rds", sep="");
    }else if(type == "refseq"){
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.path <- paste(lib.path, data.org, "/entrez_refseq.rds", sep="");
    }else if(type == "emblgene"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_gene.rds", sep="");
    }else if(type == "tair"){
      q.mat <- do.call(rbind, strsplit(q.vec, "\\."));
      q.vec <- q.mat[,1];
      db.path <- paste(lib.path, data.org, "/tair.rds", sep="");
    }else if(type == "embltranscript"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_transcript.rds", sep="");
    }else if(type == "emblprotein"){
      db.path <- paste(lib.path, data.org, "/entrez_embl_protein.rds", sep="");
    }else if(type == "orfid"){ # only for yeast
      db.path <- paste(lib.path, data.org, "/entrez_orf.rds", sep="");
    }else if(type == "flybase"){
      db.path <- paste(lib.path, data.org, "/entrez_flybase.rds", sep="");
    }else if(type == "string"){ 
      db.path <- paste(lib.path, data.org, "/entrez_string.rds", sep="")
    }else if(type == "ecogene"){ # only for ecoli
      db.path <- paste(lib.path, data.org, "/entrez_ecogene.rds", sep="")
    }else if(type == "uniprot"){
      db.path <- paste(lib.path, data.org, "/entrez_uniprot.rds", sep="")
    }else if(type == "paelocus"){
      db.path <- paste(lib.path, data.org, "/entrez.rds", sep="")
    }else{
      print("Unknown data type");
      return(0);
    }
    db.map <-  readRDS(db.path);
    hit.inx <- match(q.vec, db.map[, "accession"]);
    entrezs <- db.map[hit.inx, ];
  }
  entrezs = entrezs[c(2,1)];
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
    GeneAnotDB <-doProteinIDMapping(rownames(gene.mat), type);
    
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
