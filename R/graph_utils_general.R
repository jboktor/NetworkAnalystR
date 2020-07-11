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
    net.stats[i,] <- c(vcount(g),ecount(g),sum(V(g)$name %in% seed.proteins));
    #print(net.stats[i, ]);
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
  }else if(algo == "fr"){
    pos.xy <- layout_with_fr(g, dim=3, niter=500)
    p <- stack(data.frame(pos.xy))
    p[,1] = rescale2NewRange(p[,1], -1, 1)
    pos.xy <- unstack(p)
  }else if(algo == "kk"){
    pos.xy <- layout_with_kk(g, dim=3, maxiter=500)
    p <- stack(data.frame(pos.xy))
    p[,1] = rescale2NewRange(p[,1], -1, 1)
    pos.xy <- unstack(p)
  }
  pos.xy;
}

UpdateNetworkLayout <- function(algo, filenm){
  current.net <- ppi.comps[[current.net.nm]];
  #pos.xy <- PerformLayOut_mem(current.net.nm, algo);
  pos.xy <- PerformLayOut(current.net.nm, algo);
  nms <- V(current.net)$name;
  nodes <- vector(mode="list");
  if(algo %in% c("fr", "kk")){
    for(i in 1:length(nms)){
      nodes[[i]] <- list(
        id=nms[i],
        x=pos.xy[i,1],
        y=pos.xy[i,2],
        z=pos.xy[i,3]
      );
    }
  }else{
    for(i in 1:length(nms)){
      nodes[[i]] <- list(
        id=nms[i], 
        x=pos.xy[i,1], 
        y=pos.xy[i,2]
      );
    }
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
    if(netUploadU == 1){
      qnums <- psize;
    }else{
      hits <- seed.proteins %in% path.ids;
      qnums <- sum(hits);
    }
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
  pvall <<- pval.vec
  if(length(pval.vec)>1){
    ord.inx <- order(pval.vec, decreasing=F);
    community.vec <- community.vec[ord.inx];
    qnum.vec <- qnum.vec[ord.inx];
    ord.inx <- order(qnum.vec, decreasing=T);
    community.vec <- community.vec[ord.inx];
  }
  
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
  current.net.nm <<- net.nm
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
  
  if(length(V(g)$name)>100){
    modules = FindCommunities("walktrap", FALSE);
  }else{
    modules = "NA"
  }
  
  # covert to json
  require(RJSONIO);
  netData <- list(nodes=nodes, edges=edge.mat, org=data.org, analType=anal.type, naviString = "network", modules=modules);
  partialToBeSaved <<- c(partialToBeSaved, c(filenm))
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
