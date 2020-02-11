#' @title ReadTabExpressData
#' @description read tab delimited file
#' @param fileName file location
#' @return return a list (data.name, data.frame, meta.data)
#' @details Can have many classes, stored in meta.info (starting with #)
#' @rdname ReadTabExpressData
#' @export 
ReadTabExpressData <- function(fileName) {
  
  dataSet <- .readTabData(fileName);
  
  # rename data to data.orig
  int.mat <- dataSet$data;
  dataSet$data <- NULL;
  
  msg <- paste("a total of ", ncol(int.mat), " samples and ", nrow(int.mat), " features were found");
  
  # remove NA, null
  row.nas <- apply(is.na(int.mat)|is.null(int.mat), 1, sum);
  good.inx<- row.nas/ncol(int.mat) < 0.5;
  if(sum(!good.inx) > 0){
    int.mat <- int.mat[good.inx,];
    msg <- c(msg, paste("removed ", sum(!good.inx), " features with over 50% missing values"));
  }
  # remove constant values
  filter.val <- apply(int.mat, 1, IQR, na.rm=T);
  good.inx2 <- filter.val > 0;
  if(sum(!good.inx2) > 0){
    int.mat <- int.mat[good.inx2,];
    msg <- c(msg, paste("removed ", sum(!good.inx2), " features with constant values"));
  }
  
  if(nrow(int.mat) > 5000){
    filter.val <- filter.val[good.inx2];
    rk <- rank(-filter.val, ties.method='random');
    
    var.num <- nrow(int.mat);
    kept.num <- 0.95*var.num;
    int.mat <- int.mat[rk < kept.num, ];
    # msg <- c(msg, paste("removed 5% features with near-constant values"));
  }
  
  minVal <- min(int.mat, na.rm=T);
  na.inx <- is.na(int.mat);
  if(sum(na.inx) > 0){
    int.mat[na.inx] <- minVal/2;
    # msg <- c(msg, "the remaining", sum(na.inx), "missing variables were replaced with data min");
  }
  current.msg <<- paste(msg, collapse="; ");
  data.proc <- RemoveDuplicates(int.mat, "mean", quiet=T);
  dataSet$smpl.num <- ncol(data.proc);
  
  # save processed data for download user option
  write.csv(data.proc, file="data_processed.csv");
  saveRDS(data.proc, "data.proc.rds");
  
  dataSet <<- dataSet;
  return (1);
}

#' @title .readDataTable
#' @description read tab delimited file
#' @param fileName file location
#' @return return a list (data.name, data.frame, meta.data)
#' @details Can have many classes, stored in meta.info (starting with #)
#' @rdname .readDataTable
#' @export 
.readDataTable <- function(fileName){
  if(length(BiocGenerics::grep('\\.zip$',fileName,perl=TRUE))>0){
    fileName <- utils::unzip(fileName);
    if(length(fileName) > 1){
      # test if "__MACOSX" or ".DS_Store"
      osInx <- BiocGenerics::grep('MACOSX',fileName,perl=TRUE);
      if(length(osInx) > 0){
        fileName <- fileName[-osInx];
      }
      dsInx <- BiocGenerics::grep('DS_Store',fileName,perl=TRUE);
      if(length(dsInx) > 0){
        fileName <- fileName[-dsInx];
      }
      dat.inx <- BiocGenerics::grep(".[Tt][Xx][Tt]$", fileName);
      if(length(dat.inx) != 1){
        current.msg <<- "More than one text files (.txt) found in the zip file.";
        return(0);
      }
    }
  }
  dat <- try(data.table::fread(fileName, header=TRUE, check.names=FALSE, data.table=FALSE));
  if(class(dat) == "try-error"){
    # try to use "tr" to remove double return characters
    trFileName <- BiocGenerics::paste("tr -d \'\\r\' <", fileName);
    dat <- try(data.table::fread(trFileName, header=TRUE, check.names=FALSE, data.table=FALSE));
    if(class(dat) == "try-error"){
      ipred::print("Using slower file reader ...");
      formatStr <- substr(fileName, nchar(fileName)-2, nchar(fileName))
      if(formatStr == "txt"){
        dat <-try(readr::read.table(fileName,header=TRUE,comment.char = "", check.names=F, as.is=T));
      }else{ # note, read.csv is more than read.table with sep=","
        dat <-try(readr::read.csv(fileName,header=TRUE,comment.char = "", check.names=F, as.is=T));
      }  
    }
  }
  return(dat);
}

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
