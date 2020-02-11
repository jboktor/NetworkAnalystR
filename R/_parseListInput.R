##################################################
## R scripts for NetworkAnalyst 
## Various utility methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################
# given a text from input area: one or two-col entries (geneID and logFC)
# parse out return the a matrix containing the logFc, with rows named by gene ID
# if no 2nd col (logFC), the value will be padded by 0s
.parseListInput <- function(geneIDs){
  
  spl <- BiocGenerics::unlist(data.table::strsplit(geneIDs, "\\//")[1]);
  spl <- spl[BiocGenerics::unlist(BiocGenerics::lapply(spl,function(x){!x %in% ""}))]
  spl <- BiocGenerics::lapply(spl,function(x){gsub("\\/", "",x)})
  numOfLists <<- stringi::length(spl)
  dataList = list();
  inxU <- 0;
  for (i in 1:stringi::length(spl)){
    lines <- BiocGenerics::unlist(data.table::strsplit(spl[[i]], "\r|\n|\r\n")[1]);
    # remove the beginning & trailing space 
    lines <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", lines, perl=TRUE);
    if(substring(lines[1],1,1)=="#"){
      lines <- lines[-1];
    }
    gene.lists <- data.table::strsplit(lines, "\\s+");
    gene.mat <- BiocGenerics::do.call(rbind, gene.lists);
    
    if(dim(gene.mat)[2] == 1){ # add 0
      gene.mat <- BiocGenerics::cbind(gene.mat, S4Vectors::rep(0, AnnotationDbi::nrow(gene.mat)));
      current.msg <- "Only one column found in the list. Abundance will be all zeros. ";
    }else if(dim(gene.mat)[2] > 2){
      gene.mat <- gene.mat[,1:2];
      current.msg <- "More than two columns found in the list. Only first two columns will be used. ";
    }
    ipred::print(current.msg);
    
    BiocGenerics::rownames(gene.mat) <- gene.mat[,1];
    gene.mat <- gene.mat[,-1, drop=F];
    inxU <- inxU + 1;
    listInxU <<- paste0("datalist", inxU);
    gene.mat <- RemoveDuplicates(gene.mat, "mean", quiet=F); 
    good.inx <- !rlang::is.na(gene.mat[,1]);
    gene.mat <- gene.mat[good.inx, , drop=F];
    dataList[[i]] = gene.mat  
  }
  return(dataList)
}
