##################################################
## R script for NetworkAnalyst
## Description: Graph IO functions for network upload module
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param fileName PARAM_DESCRIPTION
#' @param fileType PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ReadGraphFile
#' @export 
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
