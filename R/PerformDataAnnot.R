# read in the data and perform
# gene ID mapping using built in libraries
# matchMin is minimal matched probe (%)
# idType: INVEX supported ID types
# lvlOpt: "NA" to keep original, other values will merge original ID to entrez gene IDs
# return the total matched gene number
# note: unmapped IDs will be retained as 
# original label (i.e. intergenic regions) in further analysis
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param org PARAM_DESCRIPTION
#' @param dataType PARAM_DESCRIPTION
#' @param idType PARAM_DESCRIPTION
#' @param lvlOpt PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformDataAnnot
#' @export 
PerformDataAnnot <- function(org, dataType, idType, lvlOpt){
  data.org <<- org;
  SetInitLib(org)
  
  dataSet$type <- dataType;
  dataSet$id.orig <- dataSet$id.current <- idType;
  dataSet$annotated <- F;
  # should not contain duplicates, however sanity check
  data.proc <- readRDS("data.proc.rds");
  dataSet$data.anot <- data.proc;
  
  if (org != 'NA' & idType != 'NA'){
    feature.vec <- rownames(data.proc);
    anot.id <- doAnnotation(feature.vec, idType);
    
    #dataSet$annotation <- anot.id; 
    saveRDS(anot.id, "annotation.rds");
    
    hit.inx <- !is.na(anot.id);
    matched.len <- sum(hit.inx);
    perct <- round(matched.len/length(feature.vec),3)*100;
    thresh <- 0.1 # previous value of 0.25 is causing challenges 
    #for datasets like Ppromelas with low annotation quality
    if (matched.len < length(feature.vec)*thresh){
      current.msg <<- paste('Only ', perct, '% ID were matched. You may want to choose another ID type or use default.', sep=""); 
    } else {
      current.msg <<- paste("ID annotation: ", "Total [", length(anot.id), 
                            "] Matched [", matched.len, "] Unmatched [", sum(!hit.inx),"]", collapse="\n");    
      
      if (lvlOpt != 'NA' | idType == "entrez"){
        # do actual summarization to gene level
        
        matched.entrez <- anot.id[hit.inx];
        data.anot <- data.proc[hit.inx,];
        rownames(data.anot) <- matched.entrez;
        current.msg <<- paste(current.msg, "Data is now transformed to gene-level (Entrez) expression.");
        
        dataSet$data.anot <- RemoveDuplicates(data.anot, lvlOpt, quiet=F);
        dataSet$id.current <- "entrez";
        dataSet$annotated <- T; 
      } else {
        current.msg <<- paste(current.msg, "No gene level summarization was performed.");
      }
    }
  } else { # no conversion will be performed
    feature.vec <- rownames(data.proc);
    anot.id = feature.vec
    perct <- 100;
    hit.inx <- !is.na(anot.id);
    matched.len <- length(feature.vec); # dummies
    minLvl <- 1;
    current.msg <<- paste("No annotation was performed. Make sure organism and gene ID are specified correctly!"); 
  }
  # need to save the ids (mixed gene annotation and original id) 
  # in case, users needs to keep unannotated features
  # this need to be updated to gether with data from now on
  dataSet$data.norm <- dataSet$data.anot;
  dataSet <<- dataSet;
  
  saveRDS(dataSet$data.anot, file="orig.data.anot"); # keep original copy, not in mem
  
  totalCount =  sum(colSums(dataSet$data.anot));
  avgCount = sum(colSums(dataSet$data.anot))/ ncol(dataSet$data.anot);
  minCount = min(colSums(dataSet$data.anot))
  maxCount = max(colSums(dataSet$data.anot))
  
  if(length(dataSet$meta.info)==1){
    lvls = paste(levels(dataSet$meta.info[,1]),collapse="; ")
  }else{
    conc1 = paste0("<b>", colnames(dataSet$meta.info)[1], "</b>", ": ", paste(levels(dataSet$meta.info[,1]), collapse="; "))
    conc2 = paste0("<b>", colnames(dataSet$meta.info)[2], "</b>", ": ", paste(levels(dataSet$meta.info[,2]), collapse="; "))
    lvls = paste("Two factors found -", conc1, conc2)
  }
  summaryVec <<- c(matched.len, perct, length(anot.id), sum(!hit.inx), ncol(dataSet$data.anot), ncol(dataSet$meta.info), sprintf("%4.2e", signif(totalCount ,3)), sprintf("%4.2e",signif(avgCount, 3)), sprintf("%4.2e",signif(minCount, 3)), sprintf("%4.2e",signif(maxCount,3)), lvls)  
  return(matched.len);   
}
