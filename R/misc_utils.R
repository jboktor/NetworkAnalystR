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
  return(length(dataSet$sig.mat))}


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
  return (colors)}

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
  return (query)}

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
  return (query)}

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
  ppi.res <- dbSendQuery(ppi.db, statement);
  ppi.res <- fetch(ppi.res, n=-1); # get all records
  dbDisconnect(ppi.db); cleanMem();

  # remove dupliated edges
  ppi.res <- ppi.res[!duplicated(ppi.res$row_id),]
  return(ppi.res)}

cleanMem <- function(n=10) { for (i in 1:n) gc() }

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
  return(currentSet)}

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
  return(res)}

PrepareEnrichNet<-function(netNm, type, overlapType){
  hits <-  enr.mat[,"Hits"];
  pvals <- enr.mat[,"P.Value"];
  require(igraph);
  require(reshape);
  pvalue <- pvals;
  id <- names(pvalue);
  readRDS("current_geneset.rds") ;
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
  unlist(sapply(V(g)$name, function(x) which(x == id))) ;

  V(g)$color <- ComputeColorGradient(-log(normalize(pvalue) + min(pvalue/2)), "black", F);
  V(g)$colorw <- ComputeColorGradient(-log(normalize(pvalue) + min(pvalue/2)), "white", F);

  cnt <- hits;
  names(cnt) <- id;
  cnt2 <- cnt[V(g)$name];

  V(g)$size <- rescale(log(cnt2, base=10), 8, 32);

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
  unlist(sapply(V(bg)$name, function(x) which(x == id))) ;
  color_scale("red", "#E5C494") ;

  V(bg)$color[V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(-log(pvalue), "black", F);
  V(bg)$colorw[V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(-log(pvalue), "white", F);
  if(anal.type == "onedata"){
    tbl = dataSet$resTable
    tbl = tbl[which(doEntrez2SymbolMapping(rownames(tbl)) %in% V(bg)$name),]
    expr.val = tbl[,selectedFactorInx];
    expvals = expr.val;
    V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(expr.val, "black", T);
    V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(expr.val, "black", T);
  }else if(anal.type == "genelist" && sum(all.prot.mat[,1]) != 0){
    tbl = all.prot.mat
    gene.nms = V(bg)$name[which(!V(bg)$name %in% rownames(enr.mat))]
    tbl = tbl[which(rownames(tbl) %in% gene.nms),]
    expr.val = as.vector(tbl);
    expvals = expr.val;
    V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(expr.val, "black", T);
    V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(expr.val, "black", T);

  }else if(anal.type =="metadata"){
    tbl = meta.mat.all
    tbl = as.matrix(tbl[which(doEntrez2SymbolMapping(rownames(tbl)) %in% V(bg)$name),])
    expr.val = tbl[,1];
    expvals = expr.val;
    V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(expr.val, "black", T);
    V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- ComputeColorGradient(expr.val, "black", T);
  }else{
    expvals <- rep(0,length(V(bg)$color));
    V(bg)$color[!V(bg)$name %in% rownames(enr.mat)] <- "#00FFFF";
    V(bg)$colorw[!V(bg)$name %in% rownames(enr.mat)] <- "#668B8B"
  }

  V(bg)$size <- rescale(log(cnt2, base=10), 8, 24);

  # layout
  pos.xy <- layout.auto(bg);

  # now create the json object
  bnodes <- vector(mode="list");
  node.nms <- V(bg)$name;
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
      exp= expvals[i],
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
      0;
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
    0;
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
  return(res)}

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
  return(com.symbols)}

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
  return(com.symbols)}

rescale <- function(x, from, to){
  (x - min(x)) / max(x - min(x)) * (to - from) + from
}

GetDataListNames <- function(){
  return(names(mdata.all))}

normalize <- function(x)
{ 
  return((x- min(x)) /(max(x)-min(x)))
}
