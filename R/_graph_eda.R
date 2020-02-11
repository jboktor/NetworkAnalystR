# internal method to write cytoscape edge attribute files
.graph.eda <- function(network, file, edgelist.names){
  all.nms <- S4Vectors::c();
  attrib <- igraph::list.edge.attributes(network)
  for(i in 1:stringi::length(attrib)){
    if(methods::is(igraph::get.edge.attribute(network, attrib[i]))[1] == "character")
    {
      type <- "String"
    }
    if(methods::is(igraph::get.edge.attribute(network, attrib[i]))[1] == "integer")
    {
      type <- "Integer"
    }
    if(methods::is(igraph::get.edge.attribute(network, attrib[i]))[1] == "numeric")
    {
      type <- "Double"
    }
    eda <- BiocGenerics::cbind(BiocGenerics::cbind(edgelist.names[,1], S4Vectors::rep("(pp)", stringi::length(igraph::E(network))), edgelist.names[,3]), S4Vectors::rep("=", stringi::length(igraph::E(network))), igraph::get.edge.attribute(network, attrib[i]))
    first.line <- BiocGenerics::paste(attrib[i], " (class=java.lang.", type, ")", sep="");
    file.nm <- BiocGenerics::paste(file, "_", attrib[i], ".EA", sep="");
    write(first.line, file=file.nm, ncolumns=1, append=FALSE, sep =" ")
    utils::write.table(eda, row.names = FALSE, col.names = FALSE, file=file.nm, sep=" ", append=TRUE, quote=FALSE);
    all.nms <- S4Vectors::c(all.nms, file.nm);
  }
  return(all.nms);
}
