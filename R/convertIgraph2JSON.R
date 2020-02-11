#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param net.nm PARAM_DESCRIPTION
#' @param filenm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname convertIgraph2JSON
#' @export 
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
