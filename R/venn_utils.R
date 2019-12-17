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
  list() ;
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
  return(1)}

PrepareSelVennData<-function(selectedNms){
  list() ;
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
    Prepare2Venn(sel.dats) ;
  }else if(length(sel.dats) == 3){
    Prepare3Venn(sel.dats) ;
  }else if(length(sel.dats) == 4){
    Prepare4Venn(sel.dats) ;
  }
  venn.list.up <<- sel.dats;
  venn.genenb.up <<- venn.genenb
  return(1)}


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
  return(length(venn.list))}

GetSelectedDataNames <- function(){
  return(paste(names(venn.list), collapse=";"))}

GetSelectedDataGeneNumber<- function(){
  return(paste(venn.genenb, collapse=";"))}

GetSelectedDataNumberUpdated <- function(){
  return(length(venn.list.up))}

GetSelectedDataNamesUpdated <- function(){
  return(paste(names(venn.list.up), collapse=";"))}

GetSelectedDataGeneNumberUpdated<- function(){
  return(paste(venn.genenb.up, collapse=";"))}


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
  return(paste(unique(sym.vec), collapse="||"))}

PerformVennEnrichment <- function(file.nm, fun.type){
  res <- PerformEnrichAnalysis(file.nm, fun.type, venn.genes);
  return(res)}
