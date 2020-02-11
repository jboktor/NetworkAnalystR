.prepareSigProteinJSON <- function(){
  if(anal.type == "genelist"){
    result.list <- .prepareListSeeds();
  }else{ # single expression data or meta.mat
    result.list <- .prepareExpressSeeds();
  }
  return(result.list);
}
