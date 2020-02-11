#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param imgNm PARAM_DESCRIPTION
#' @param dpi PARAM_DESCRIPTION, Default: 72
#' @param format PARAM_DESCRIPTION, Default: 'png'
#' @param factor PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PlotLibSizeView
#' @export 
PlotLibSizeView<-function(imgNm,dpi=72, format="png",factor){
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  dpi=as.numeric(dpi)
  data_bef<-data.matrix(dataSet$data.anot);
  
  smpl.sums <- colSums(data_bef);
  
  library("ggplot2")
  data_bef<-data.matrix(dataSet$data.anot);
  smpl.sums <- colSums(data_bef);
  names(smpl.sums) <- colnames(data_bef);
  sampleNms = names(smpl.sums)
  df = data.frame(count=smpl.sums,ind=colnames(data_bef))
  
  if(length(dataSet$meta.info) == 2){
    Factor1 = as.vector(dataSet$meta.info[,1])
    factor1Nm = colnames(dataSet$meta.info)[1]
    conv = data.frame(ind=sampleNms, Factor1=Factor1)
    colnames(conv) = c("ind", factor1Nm)
    df1 = merge(df, conv, by="ind")
    Factor2 = as.vector(dataSet$meta.info[,2])
    factor2Nm = colnames(dataSet$meta.info)[2]
    conv = data.frame(ind=sampleNms, Factor2=Factor2)
    colnames(conv) = c("ind", factor2Nm)
    df1 = merge(df1, conv, by="ind")
    df2 <- melt(df1, measure.vars=c(factor1Nm,factor2Nm))
    colnames(df2)[4] = "Conditions"
    if(length(df2$ind)>20){
      g = ggplot(df2, aes(x = Conditions, y = count, fill=Conditions))+
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), dotsize=0.7) + ylab("Sum") + facet_grid(. ~ variable);
    }else{
      g = ggplot(df2, aes(x = Conditions, y = count, fill=Conditions, label=ind))+
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), dotsize=0.7) + geom_text_repel(force=5) + ylab("Sum") + facet_grid(. ~ variable);
    }
    width = 12
    height = 6
  }else{
    Conditions= as.character(dataSet$meta.info[,1]);
    conv = data.frame(ind=sampleNms, Conditions=Conditions)
    df1 = merge(df, conv, by="ind")
    if(length(df1$ind)>20){
      g = ggplot(df1, aes(x = Conditions, y = count, fill=Conditions))+
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), dotsize=0.7) + xlab("Sum");
    }else{
      g = ggplot(df1, aes(x = Conditions, y = count, label=ind, fill=Conditions))+
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), dotsize=0.7) + geom_text_repel(force=5) + xlab("Sum");
    }
    width = 8
    height = 6
  }
  
  Cairo(file=imgNm, width=width, height=height, unit="in", type=format, bg="white",dpi=dpi);
  
  #names(smpl.sums) <- colnames(data_bef);
  print(g);
  dev.off(); 
}
