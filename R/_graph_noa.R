# internal method to write cytoscape node attribute files
.graph.noa <- function(network, file){
  all.nms <- S4Vectors::c();
  attrib <- igraph::list.vertex.attributes(network)
  for(i in 1:stringi::length(attrib)){
    if(methods::is(igraph::get.vertex.attribute(network, attrib[i]))[1] == "character")
    {
      type <- "String"
    }
    if(methods::is(igraph::get.vertex.attribute(network, attrib[i]))[1] == "integer")
    {
      type <- "Integer"
    }
    if(methods::is(igraph::get.vertex.attribute(network, attrib[i]))[1] == "numeric")
    {
      type <- "Double"
    }
    noa <- BiocGenerics::cbind(igraph::V(network)$name, S4Vectors::rep("=", stringi::length(igraph::V(network))), igraph::get.vertex.attribute(network, attrib[i]))
    first.line <- BiocGenerics::paste(attrib[i], " (class=java.lang.", type, ")", sep="")
    file.nm <- BiocGenerics::paste(file, "_", attrib[i], ".NA", sep="");
    write(first.line, file=file.nm, ncolumns = 1, append=FALSE, sep=" ")
    utils::write.table(noa, row.names = FALSE, col.names = FALSE, file=file.nm, sep=" ", append=TRUE, quote=FALSE);
    all.nms <- S4Vectors::c(all.nms, file.nm);
  }
  return(all.nms);
}
