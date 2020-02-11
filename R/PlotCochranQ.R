#compute Cochran’s Q to help FEM/REM
# plot Q-Q plot for estimation
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param imgNm PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PlotCochranQ
#' @export 
PlotCochranQ <- function(imgNm){
  if(!performedDE){
    PerformMetaDeAnal();
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  sel.inx <- mdata.all==1;
  sel.nms <- names(mdata.all)[sel.inx];
  nbstudies <- sum(sel.inx);
  ES<-array(dim=c(nrow(inmex.meta$data),4,nbstudies));
  cls.lvls <- levels(as.factor(inmex.meta$cls.lbl));
  
  for (i in 1:nbstudies){
    data.nm <- sel.nms[i];
    dataSet <- readRDS(data.nm);
    
    fit.obj.nm <- paste(data.nm, "fit.obj", sep=".");
    fit2i <- readRDS(fit.obj.nm);
    #fit2i <- dataSet$fit.obj;
    
    n1i=length(which(dataSet$cls==cls.lvls[1]));
    n2i=length(which(dataSet$cls==cls.lvls[2]));
    ES[,,i]=effectsize(fit2i$t,((n1i*n2i)/(n1i+n2i)),(fit2i$df.prior+fit2i$df.residual));
  }
  
  Qvals <- f.Q.NA(ES[,3,],ES[,4,]);
  
  Cairo(file=imgNm, width=400, height=400, type="png", bg="white");
  # histgram
  # hist(Qvals, breaks = 50, col = "red");
  # QQ plot
  chisqq <- qchisq(seq(0, 0.9999, 0.001), df = nbstudies-1)
  tmp <- quantile(Qvals, seq(0, 0.9999, 0.001))
  qqplot(chisqq, tmp, ylab = "Quantiles of Sample", pch = "*", 
         xlab = "Quantiles of Chi > square", main = "QQ Plot")
  lines(chisqq, chisqq, lty = "dotted", col = "red")
  
  dev.off(); 
  
  Qmean <- round(mean(Qvals),5);
  return (Qmean);
}
