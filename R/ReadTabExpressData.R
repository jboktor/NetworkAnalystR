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
