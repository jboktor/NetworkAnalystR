# given a gene id, plot its expression profile as box plot
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param gene.id PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PlotSelectedGene
#' @export 
PlotSelectedGene<-function(gene.id, type){
  
  imgName <- paste("Gene_", gene.id, ".png", sep="");
  require(lattice);
  if(anal.type == "onedata"){
    ids <- rownames(dataSet$resTable);
    inx <- which(ids == gene.id);
    symb <- dataSet$sig.genes.symbols[inx]; 
    if(type== "volcano"){
      symb = "";
    }    
    if(dataSet$comp.type == "custom"){
      Cairo(file = imgName, width=280, height=320, type="png", bg="white");
      grp.nms = dataSet$grp.nms;
      inx = dataSet$cls %in% grp.nms;
      cls = dataSet$cls[inx]
      dat = dataSet$data.norm[,inx];
      myplot <- bwplot(dat[gene.id,] ~ as.character(cls), fill="#0000ff22", scales=list(x=list(rot=30)),
                       xlab="Class", ylab="Expression Pattern", main=symb);
    }else if(length(dataSet$sec.cls)>1){
      out.fac <- as.character(dataSet$sec.cls)
      in.fac <- as.character(dataSet$fst.cls)
      xlab = colnames(dataSet$meta.info[,1]);
      
      Cairo(file = imgName, dpi=72, width=320, height=320, type="png", bg="white");
      #ylim.ext <- GetExtendRange(dataSet$data.norm[gene.id, ], 12);
      layout <- c(2, 1);
      myplot<- bwplot(dataSet$data.norm[gene.id, ] ~ in.fac | out.fac, 
                      xlab="Factors", ylab="Expression Pattern", main=symb, scales=list(x=list(rot=30)),
                      fill="#0000ff22", layout=layout);
    }else{
      Cairo(file = imgName, width=280, height=320, type="png", bg="white");
      
      myplot <- bwplot(dataSet$data.norm[gene.id,] ~ as.character(dataSet$cls), fill="#0000ff22", scales=list(x=list(rot=30)),
                       xlab="Class", ylab="Expression Pattern", main=symb);
    }
    
  }else{ # metadata
    
    inmex.meta <- readRDS("inmex_meta.rds");
    if(inmex.meta$id.type == "entrez"){
      symb <- inmex.meta$gene.symbls[gene.id];
    }else{
      symb <- gene.id;
    }
    num <- sum(mdata.all == 1);
    # calculate width based on the dateset number
    if(num == 1){
      Cairo(file = imgName, width=280, height=320, type="png", bg="white");
      myplot <- bwplot(inmex.meta$plot.data[gene.id,] ~ as.character(inmex.meta$cls.lbl), fill="#0000ff22",
                       xlab="Class", ylab="Expression Pattern", main=symb, scales=list(x=list(rot=30)))
    }else{
      # calculate layout
      if(num < 6){
        layout <- c(num, 1);
        height=320;
        width=160*num;
      }else{
        rn <- round(num/2);
        layout <- c(rn, 2);
        height=500;
        width=160*rn;
      }
      
      Cairo(file = imgName, width=width, height=height, type="png", bg="white");
      data.lbl <- as.character(inmex.meta$data.lbl);
      data.lbl <- substr(data.lbl, 0, nchar(data.lbl)-4);
      
      # get counts in each data, same order as a levels
      counts <- table(data.lbl);
      # back to factor 
      data.lbl <- factor(data.lbl);
      
      # get new lbls to cut potential long names, and add sample numbers
      nlbls <- data.lbl;
      levels(nlbls) <- abbreviate(levels(nlbls),9);
      nlbls <- paste(levels(nlbls), "( n=", as.vector(counts), ")");
      # update labels
      data.lbl <- factor(data.lbl, labels=nlbls);
      # some time the transformed plot.data can switch class label, use the original data, need to be similar scale
      myplot <- bwplot(inmex.meta$plot.data[gene.id,] ~ as.character(inmex.meta$cls.lbl) | data.lbl, 
                       xlab="Datasets", ylab="Expression Pattern", main=symb, scales=list(x=list(rot=30)),
                       fill="#0000ff22", layout=layout);
    }
  }
  print(myplot); 
  dev.off();
}
