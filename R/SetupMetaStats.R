# for gene-level meta-analysis
# function to set up results combining individual data analysis
# as well as to prepare for GO analysis
# no return, as set global 
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param BHth PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname SetupMetaStats
#' @export 
SetupMetaStats <- function(BHth){
  
  GlobalCutOff$BHth <<- BHth;
  #all common genes
  inmex.meta <- readRDS("inmex_meta.rds");
  gene.ids <- rownames(inmex.meta$data);
  # meta.sig genes
  metade.genes <- rownames(meta.mat);
  
  # setup individual sig genes & stats
  # that overlap with meta.sig
  inmex.de <- list();
  
  pval.mat <- fc.mat <- matrix(nrow=nrow(meta.mat), ncol=sum(mdata.all==1));
  
  for(i in 1:length(inmex.ind)){
    de.res <- inmex.ind[[i]];
    
    hit.inx <- de.res[,2] <= BHth;
    hit.inx <- which(hit.inx); # need to get around NA
    inmex.de[[i]] <- rownames(de.res)[hit.inx];
    
    # only choose the genes that are also meta sig genes from in
    # individual analysis for display
    de.res <- de.res[metade.genes,];
    
    fc.mat[,i] <- de.res[,1];
    pval.mat[,i] <- de.res[,2];
  }
  names(inmex.de) <- names(inmex.ind);
  
  # calculate gain/loss
  deindst <- unique(unlist(inmex.de));
  gains=metade.genes[which(!(metade.genes %in% deindst))];
  losses=deindst[which(!(deindst %in% metade.genes))];
  all.de <- cbind(gene.ids %in% metade.genes, gene.ids %in% deindst);
  colnames(all.de) <- c("Meta-DE", "Individual-DE");
  vennC <- getVennCounts(all.de);
  if(inmex.meta$id.type == "entrez"){ 
    names(metade.genes) <- inmex.meta$gene.symbls[metade.genes];
    names(gains) <- inmex.meta$gene.symbls[gains];
    names(losses) <- inmex.meta$gene.symbls[losses];
  }
  
  # de genes from individual 
  de.len <- sapply(inmex.de, length);
  stat <- c(length(metade.genes), de.len);
  names(stat) <- c("Meta", substr(names(inmex.de), 0, nchar(names(inmex.de))-4));
  meta.stat <- list(
    stat = stat,
    de = metade.genes,
    idd = gains,
    loss = losses,
    venn = vennC
  );
  
  fc.mat <<- fc.mat;
  pval.mat <<- pval.mat;
  inmex.de <<- inmex.de;
  meta.stat <<- meta.stat;
  
  # save the result
  if(inmex.meta$id.type == "entrez"){ # row name gene symbols
    metade.nms <- inmex.meta$gene.symbls[metade.genes];
    res <- cbind(EntrezID=metade.genes, Name=metade.nms, meta.mat);
  }else{
    res <- cbind(ID=metade.genes, meta.mat);
  }
  write.csv(res, file=paste("meta_sig_genes_", inmex.method, ".csv", sep=""), row.names=F);
}
