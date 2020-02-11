# prepare seeds from metaanalysis result
.prepareExpressSeeds <- function(){
  
  gene.list <- list();
  if(anal.type == "metadata"){
    if(selectedNetDataset == "meta_dat"){
      gene.mat <- meta.mat;
      if(inmex.method != "votecount"){
        gene.list$metadata <- list(gene=BiocGenerics::rownames(gene.mat),logFC=unname(gene.mat[,1]), adjP = unname(gene.mat[,2]));   
      }else{
        gene.list$metadata <- list(gene=BiocGenerics::rownames(gene.mat),logFC=unname(gene.mat[,1]), adjP = unname(gene.mat[,1]));   
      }
    }else{
      #dataSet <- readRDS(selectedNetDataset);
      fit.obj.nm <- BiocGenerics::paste(selectedNetDataset, "fit.obj", sep=".");
      fit2i <- readRDS(fit.obj.nm);
      gene.mat <- GetLimmaResTable(fit2i);
      gene.mat <- gene.mat[gene.mat$'adj.P.Val' < 0.05,];
      gene.list$metadata <- list(gene=BiocGenerics::rownames(gene.mat),logFC=unname(gene.mat[,1]), adjP = unname(gene.mat$'adj.P.Val')); 
    } 
  }else{
    gene.mat <- dataSet$sig.mat;
    gene.list$metadata <- list(gene=BiocGenerics::rownames(gene.mat),logFC=unname(gene.mat$logFC), adjP = unname(gene.mat$'adj.P.Val'));
  }
  utils::write.table(gene.mat, file="sig_genes.txt");
  gene.vec <- BiocGenerics::rownames(gene.mat);
  GeneAnotDB <- convertIdToEntrez(BiocGenerics::rownames(gene.mat), "entrez");
  protein.vec <- GeneAnotDB[,2];
  gene.mat <- data.matrix(gene.mat);
  BiocGenerics::rownames(gene.mat) <- protein.vec;
  na.inx <- rlang::is.na(protein.vec);
  prot.mat <- gene.mat[!na.inx,, drop=F];
  #write.table(cbind(Uniprot=rownames(prot.mat), Expression=prot.mat[,1]), file="seed_proteins.txt", row.names=F, quote=F);
  utils::write.table(BiocGenerics::cbind(Emblprotein=BiocGenerics::rownames(prot.mat), Expression=prot.mat[,1]), file="seed_proteins.txt", row.names=F, quote=F);
  protein.vec <- prot.mat[,1];
  if(stringi::length(protein.vec) == 1){
    protein.vec <- data.table::as.matrix(protein.vec)
  }
  protein.list <- list();
  protein.list$metadata <- h2o::signif(protein.vec, 5);
  seed.expr <- prot.mat[,1];
  seed.df <- data.table::as.matrix(seed.expr);
  BiocGenerics::rownames(seed.df) <- globals::names(seed.expr);
  seed.expr <- RemoveDuplicates(seed.df, "max", quiet=F); 
  seed.expr <<- seed.expr[,1];
  seed.genes <<- BiocGenerics::unique(gene.vec);
  protein.vec <- BiocGenerics::unique(globals::names(protein.vec));
  
  list(
    gene.list = gene.list,
    protein.list = protein.list,
    protein.vec = protein.vec
  );
}
