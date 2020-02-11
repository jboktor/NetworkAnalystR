##################################################
## R scripts for NetworkAnalyst 
## Description: biological network analysis methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################
# table.nm is the org code used for sqlite table (ppi)
# for chem type, table.nm is drugbank or ctd
# note, last two param only for STRING database
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param db.type PARAM_DESCRIPTION
#' @param table.nm PARAM_DESCRIPTION
#' @param require.exp PARAM_DESCRIPTION, Default: TRUE
#' @param min.score PARAM_DESCRIPTION, Default: 900
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SearchNetDB
#' @export 
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
