# given dataSet Name, sample name, and class name, do update
# note, for multiple #class, this set which one to use in the subsequent steps
# last one wins
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataName PARAM_DESCRIPTION
#' @param clsLbl PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname UpdateSampleInfo
#' @export 
UpdateSampleInfo<-function(dataName, clsLbl){
  
  print("updating sample info .... ");
  
  if(!exists("class.vec")){
    print("Could not find class label list!");
    return(0);
  }
  
  if(!exists("smpl.vec")){
    print("Could not find sample name list!");
    return(0);
  }
  
  if(length(class.vec) < 2){
    current.msg <<- "Add least two groups required!";
    return(0);
  }
  
  if(sum(class.vec != 'NA') < 2){
    current.msg <<- "Cannot be less than 2 groups";
    return(0);
  }
  
  if(dataSet$name != dataName){
    dataSet <- readRDS(dataName);
  }
  org.lvl.len <- length(levels(dataSet$meta.info[[clsLbl]]));
  if(org.lvl.len < length(class.vec)){
    current.msg <<- "You can not add new groups";
    return(0);
  }else if(org.lvl.len > length(class.vec)){
    current.msg <<- "To exclude a group, replace it with NA.";
    return(0);
  }
  
  # first update the meta info
  cls <- dataSet$meta.info[[clsLbl]];
  levels(cls) <- class.vec;
  
  data <- dataSet$data.orig;
  meta.info <- dataSet$meta.info;
  
  if(any(levels(cls) == 'NA')){
    rt.inx <- cls != 'NA';
    data <- data[,rt.inx];
    
    # also update the whole meta-info
    meta.info <- meta.info[rt.inx,,drop=FALSE];
    cls <- cls[rt.inx];
  }
  
  # need to re-construct the class, so that the level order  
  # are always alphabetic
  meta.info[[clsLbl]] <- factor(as.character(cls));
  
  # note, sample names could be removed (together with cls) as the whole row
  hit.inx <- colnames(data)%in%smpl.vec;
  dataSet$data.orig <- data[,hit.inx];
  
  # make sure the factor levels also dropped
  for(i in 1:length(meta.info)){
    meta.info[[i]] <- factor(meta.info[[i]][hit.inx]);
  }
  
  dataSet$meta.info <- meta.info;
  dataSet$cls <-  dataSet$meta.info[[clsLbl]];
  RegisterData(dataSet);
  gc();
  return(1);
}
