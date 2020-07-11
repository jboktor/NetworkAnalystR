##################################################
## R script for NetworkAnalyst
## Description: Graph IO functions for network upload module
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################

ReadGraphFile <- function(fileName, fileType, org, idType) {
  data.org<<- org
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
  current.msg <<- "Successfully parsed your graph file!";
  if(idType == "NA"){
    nms <- V(graphX)$name;
    if(length(nms)<1){
      nms <- V(graphX)$id;
      graphX = set_vertex_attr(graphX, "name", value=nms)
    }
    node.data = data.frame(nms, nms);
    seed.proteins <<- nms;
  }else{
    nms <- V(graphX)$name;
    if(length(nms)<1){
      nms <- V(graphX)$id;
      graphX = set_vertex_attr(graphX, "name", value=nms)
    }
    entrezs = convertIdToEntrez(nms, idType)
    nms = entrezs[,"gene_id"]
    symbols = doEntrez2SymbolMapping(nms)
    node.data = data.frame(nms, symbols);
    seed.proteins <<- entrezs;
    V(graphX)$name = nms;
  }
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
  data.idType <<- idType; 
  convertIgraph2JSONFromFile(net.nm, "networkanalyst_0.json", data.idType);
  return(1);
}


# create igraph from the edgelist saved from graph DB
# and decompose into subnets

convertIgraph2JSONFromFile <- function(net.nm, filenm, idType){
  
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
      label=ppi.net$node.data[i,2],
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
  
  netUploadU <<-1
  if(length(V(g)$name)>100){
    modules = FindCommunities("walktrap", FALSE);
  }else{
    modules = "NA"
  }
  
  # covert to json
  require(RJSONIO);
  netData <- list(nodes=nodes, edges=edge.mat, idType=idType, analType=anal.type, org=data.org, naviString = "network", modules=modules);
  dataSet$jsonNms$network <<- filenm
  if(!filenm %in% partialToBeSaved){
    partialToBeSaved <<- c(partialToBeSaved, c(filenm))
  }
  
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
