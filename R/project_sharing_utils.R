PrepareSignatureOfNetworkAnalyst <- function(){
  
  if(anal.type == "genelist"){
    signature.gene <- dataSet$sig.mat;
    
    signature.gene.org <- data.org;
    save(signature.gene, signature.gene.org, file="RShare_networkanalyst.RData");  
  }else if(anal.type == "onedata"){
    if(!file.exists("ExpressResT.rda")){
      return("-1");
    }
    resT <- readRDS("ExpressResT.rda");
    if(exists("P.Value", where=resT)){
      signature.gene <- as.matrix(resT$P.Value);
    }else if(exists("PValue", where=resT)){
      signature.gene <- as.matrix(resT$PValue);
    }
    rownames(signature.gene) <- rownames(resT);
    
    signature.gene.org <- data.org;
    save(signature.gene, signature.gene.org, file="RShare_networkanalyst.RData");  
  }
  
  return(1);
}

