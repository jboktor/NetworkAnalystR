#' @title .readTabData
#' @description read tab delimited file
#' @param dataName file location
#' @return return a list (data.name, data.frame, meta.data)
#' @details read tab delimited file. Can have many classes, stored in meta.info (starting with #)
#' @rdname .readTabData
#' @export 
.readTabData <- function(dataName) {
  if(length(BiocGenerics::grep('\\.zip$',dataName,perl=TRUE))>0){
    dataName <- utils::unzip(dataName);
    if(length(dataName) > 1){
      # test if "__MACOSX" or ".DS_Store"
      osInx <- BiocGenerics::grep('MACOSX',dataName,perl=TRUE);
      if(length(osInx) > 0){
        dataName <- dataName[-osInx];
      }
      dsInx <- BiocGenerics::grep('DS_Store',dataName,perl=TRUE);
      if(length(dsInx) > 0){
        dataName <- dataName[-dsInx];
      }
      dat.inx <- BiocGenerics::grep(".[Tt][Xx][Tt]$", dataName);
      if(length(dat.inx) != 1){
        current.msg <<- "More than one text files (.txt) found in the zip file.";
        return(0);
      }
    }
  }
  
  msg <- NULL;
  # using the powerful fread function, 10 times faster, note: default return data.table, turn off
  dat1 <- .readDataTable(dataName);
  
  # look for #CLASS, could have more than 1 class labels, store in a list
  meta.info <- list();
  cls.inx <- BiocGenerics::grep("^#CLASS", dat1[,1]);
  if(length(cls.inx) > 0){ 
    for(i in 1:length(cls.inx)){
      inx <- cls.inx[i];
      cls.nm <- substring(dat1[inx, 1],2); # discard the first char #
      if(nchar(cls.nm) > 6){
        cls.nm <- substring(cls.nm, 7); # remove class
      }
      cls.lbls <- dat1[inx, -1];
      # test NA
      na.inx <- rlang::is_na(cls.lbls);
      cls.lbls[na.inx] <- "NA";
      cls.lbls <- ClearFactorStrings(cls.nm, cls.lbls);
      
      meta.info[[cls.nm]] <- cls.lbls;
    }
  }else{
    current.msg <<- "No metadata labels #CLASS found in your data!";
    return("F");
  }
  
  meta.info <- data.frame(meta.info);
  
  # now remove all comments in dat1
  # assign rownames after covert to matrix as data.frame does not allow duplicate names
  comments.inx <- BiocGenerics::grep("^#", dat1[,1]);
  dat1.nms <- dat1[-comments.inx,1];
  dat1<-dat1[-comments.inx,-1];
  dat1 <- data.matrix(dat1);
  BiocGenerics::rownames(dat1) <- dat1.nms;
  
  list(
    name= BiocGenerics::basename(dataName),
    data=dat1,
    meta.info=meta.info
  );
}
