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
  expr.vec <- abs(expr.vals)
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
  GeneAnotDB <- doProteinIDMapping(nd.nms, "entrez");
  
  entrezIDs <- GeneAnotDB[,1];
  names(entrezIDs) <- nd.nms;
  current.anot <<- entrezIDs;
  
  convertIgraph2JSON(net.nm, json.nm);
  current.net.nm <<- net.nm;
  return(1);
}

PrepareNetworkUpload <- function(net.nm, json.nm){
  
  my.ppi <- ppi.comps[[net.nm]];
  nd.nms <- V(my.ppi)$name;
  
  entrezIDs <- nd.nms
  names(entrezIDs) <- nd.nms;
  current.anot <<- entrezIDs;
  
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

GetPCSFNet <- function(){
  
  require('PCSF')
  edg = get.edgelist(overall.graph)
  edg = as.data.frame(edg)
  edg$V3 = rep(1, nrow(edg))
  colnames(edg)= c("from", "to", "cost")
  ppi = construct_interactome(edg)
  terminals = expr.vec;
  g <- PCSF(ppi, terminals, w = 5, b = 100, mu = 0.0005)
  
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
  if(vismode != "netupload"){
    prot.anots <- doProteinIDMapping(lines, id.type);
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


