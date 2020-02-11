# read in the data and perform
# gene ID mapping using built in libraries
# matchMin is minimal matched probe (%)
# return the total matched gene number
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataName PARAM_DESCRIPTION
#' @param featureType PARAM_DESCRIPTION
#' @param matchMin PARAM_DESCRIPTION, Default: 0.5
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ProcessIndData
#' @export 
ProcessIndData<-function(dataName, featureType, matchMin=0.5){
  
  dataSet <- readRDS(dataName);
  dataSet$id.type <- featureType;
  if(data.org != 'NA' & featureType != 'NA'){
    feature.vec <- rownames(dataSet$data.raw);
    minLvl <- length(feature.vec)*matchMin;
    if(featureType %in% c("entrez", "symbol", "refseq", "genbank", "emblgene", "embltranscript", "orfid", "wormbase")){
      entrez.id <- doGeneIDMapping(feature.vec, featureType);
    }else{
      entrez.id <- doProbeMapping(feature.vec, featureType);
    }   
    
    hit.inx <- which(!is.na(entrez.id));
    matched.len <- length(hit.inx);
    if(matched.len > minLvl){
      res <- "Success!";
      data.orig <- dataSet$data.raw[hit.inx,];
      matched.entrez <- entrez.id[hit.inx];
      rownames(data.orig) <- matched.entrez;
      
      # now, deal with duplicated entrez id
      # first, average duplicate rows
      ave.data <- apply(data.orig, 2, myave, matched.entrez); 
      # then removed duplicated entries
      dup.inx <- duplicated(matched.entrez);
      int.mat <- ave.data[!dup.inx,];
      # update
      dataSet$data.orig <- int.mat;
      dataSet$id.type <- "entrez";
    }
    current.msg <<- paste("ID Conversion: ", "Total [", length(entrez.id), 
                          "] Matched [", matched.len, "] Unmatched [", sum(is.na(entrez.id)),"]", collapse="\n"); 
  }else{ # no conversion will be performed
    dataSet$data.orig <- dataSet$data.raw;
    matched.len <- 9; # dummies
    minLvl <- 1;
  }
  RegisterData(dataSet);
  if(matched.len < minLvl){
    return(0);
  }
  return(matched.len);
}
