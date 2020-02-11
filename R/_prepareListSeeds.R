.prepareListSeeds <- function(){
  
  protein.list <- list();
  gene.list <- list();
  
  if(numOfLists > 1){
    if(selectedNetDataset %in% S4Vectors::c("intersect","union")){
      dataSet = list();
      dataSet$name = selectedNetDataset
      my.vec <- globals::names(mdata.all);
      com.ids <- NULL;
      list.vec <- list()
      for(i in 1:stringi::length(my.vec)){
        datSet <- readRDS(my.vec[i]);
        if(purrr::is.null(com.ids)){
          com.ids <- datSet$GeneAnotDB[,"gene_id"];
          prot.mat <- datSet$prot.mat
          list.vec[[i]] = com.ids
        }else{
          if(selectedNetDataset == "intersect"){
            com.ids <- datSet$GeneAnotDB[,"gene_id"];
            list.vec[[i]] = com.ids
            #com.ids <- intersect(com.ids, datSet$GeneAnotDB[,"gene_id"]);
          }else{
            com.ids <- BiocGenerics::union(com.ids, datSet$GeneAnotDB[,"gene_id"]);
          }
          prot.mat <- BiocGenerics::rbind(prot.mat, data.table::as.matrix(datSet$prot.mat[BiocGenerics::rownames(datSet$prot.mat) %in% com.ids,]))
        }
      }
      if(selectedNetDataset == "intersect"){
        com.ids = BiocGenerics::Reduce(intersect, list.vec)
        prot.mat <- data.table::as.matrix(datSet$prot.mat[BiocGenerics::rownames(datSet$prot.mat) %in% com.ids,])
      }else{
        com.ids <- BiocGenerics::unique(rlang::as.character(com.ids[!rlang::is.na(com.ids)])); # make sure it is interpreted as name not index
      }
      
      com.symbols <- doEntrez2SymbolMapping(com.ids);
      dataSet$GeneAnotDB = dplyr::data.frame(gene_id=com.ids, accession=com.symbols);
      dataSet$prot.mat = prot.mat;
      dataSet$sig.mat = prot.mat
      dataSet$seeds.proteins = com.ids
    }else{
      my.vec <- globals::names(mdata.all); 
      # make sure reference is the first
      inx <- BiocGenerics::which(my.vec == selectedNetDataset);
      my.vec <- my.vec[-inx];
      com.ids <- NULL;
      ids.list <- list()
      for(i in 1:stringi::length(my.vec)){
        dataSet <- readRDS(my.vec[i]);
        ids.list[[i]]=dataSet$GeneAnotDB[,"gene_id"]
      }
      dataSet <- readRDS(selectedNetDataset);
      ids <- BiocGenerics::unique(BiocGenerics::unlist(ids.list));
      com.ids <-BiocGenerics::setdiff(dataSet$GeneAnotDB[,"gene_id"], ids);
      prot.mat <- data.table::as.matrix(dataSet$prot.mat[BiocGenerics::which(BiocGenerics::rownames(dataSet$prot.mat) %in% com.ids),])
      com.symbols <- doEntrez2SymbolMapping(com.ids);
      dataSet$GeneAnotDB = dplyr::data.frame(gene_id=com.ids, accession=com.symbols);
      dataSet$prot.mat = prot.mat;
      dataSet$sig.mat = prot.mat
      dataSet$seeds.proteins = com.ids
    }
  }
  
  # return a json array object
  # each object for a single dataset its sig proteins
  meta.vec <- meta.gene.vec <- meta.seed.expr <- NULL;
  fs::file.create("seed_proteins.txt");
  GeneAnotDB <- NULL;
  
  gene.mat <- dataSet$sig.mat;
  prot.mat <- dataSet$prot.mat;
  write(BiocGenerics::paste("#DataSet:", dataSet$name),file="sig_genes.txt",append=TRUE);
  utils::write.table(dataSet$sig.mat, file="sig_genes.txt", append=TRUE);
  
  meta.gene.vec <- S4Vectors::c(meta.gene.vec, BiocGenerics::rownames(gene.mat));
  gene.list[[dataSet$name]] <- list(gene=BiocGenerics::rownames(gene.mat),logFC=unname(gene.mat[,1]));
  GeneAnotDB <- BiocGenerics::rbind(GeneAnotDB, dataSet$GeneAnotDB);
  meta.seed.expr <- S4Vectors::c(meta.seed.expr, prot.mat[,1]);
  write(BiocGenerics::paste("#DataSet:", dataSet$name),file="seed_proteins.txt",append=TRUE);
  utils::write.table(BiocGenerics::cbind(Emblprotein=BiocGenerics::rownames(prot.mat), Expression=prot.mat[,1]), file="seed_proteins.txt", row.names=F, quote=F,append=TRUE);
  protein.vec <- prot.mat[,1];
  meta.vec <- S4Vectors::c(meta.vec, globals::names(protein.vec));
  if(stringi::length(protein.vec) == 1){
    protein.vec <- data.table::as.matrix(protein.vec)
  }   
  protein.list[[dataSet$name]] <- h2o::signif(protein.vec, 3);
  
  gene.list$name <- dataSet$name;
  seed.genes <<- BiocGenerics::unique(meta.gene.vec);
  
  meta.seed.df <- data.table::as.matrix(meta.seed.expr);
  BiocGenerics::rownames(meta.seed.df) <- globals::names(meta.seed.expr);
  
  seed.expr <- RemoveDuplicates(meta.seed.df, "max", quiet=F); 
  seed.expr <<- seed.expr[,1];
  protein.vec <- BiocGenerics::unique(meta.vec);
  
  result = list(
    gene.list = gene.list,
    protein.list = protein.list,
    protein.vec = protein.vec
  );
  return(result)
}
