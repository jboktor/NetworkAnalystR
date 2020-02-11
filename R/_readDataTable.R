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
