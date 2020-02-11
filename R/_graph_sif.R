# internal function to write cytoscape .sif file
.graph.sif <- function(network, file){
  edgelist.names <- igraph::get.edgelist(network, names=TRUE)
  edgelist.names <- BiocGenerics::cbind(edgelist.names[,1], S4Vectors::rep("pp", stringi::length(igraph::E(network))), edgelist.names[,2]);
  utils::write.table(edgelist.names, row.names=FALSE, col.names=FALSE, file=BiocGenerics::paste(file, ".sif", sep=""), sep="\t", quote=FALSE)
  return(edgelist.names) 
}
