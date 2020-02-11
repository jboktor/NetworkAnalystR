#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param gene.id PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PlotSelectedGeneMeta
#' @export 
PlotSelectedGeneMeta<-function(gene.id){
  
  # first get gene symbol
  inmex.meta <- readRDS("inmex_meta.rds");
  if(inmex.meta$id.type == "entrez"){
    symb <- inmex.meta$gene.symbls[gene.id];
  }else{
    symb <- gene.id;
  }
  
  imgName <- paste("Gene_", gene.id, ".png", sep="");
  require(lattice);
  
  num <- sum(mdata.all == 1);
  # calculate width based on the dateset number
  if(num == 1){
    Cairo(file = imgName, width=280, height=320, type="png", bg="white");
    myplot <- bwplot(inmex.meta$plot.data[gene.id,] ~ as.character(inmex.meta$cls.lbl), fill="#0000ff22",
                     xlab="Class", ylab="Expression Pattern", main=symb, scales=list(x=list(rot=30)))
  }else{
    # this is a single long list 
    layout <- c(1, num);
    height <- 200*num;
    width <- 280;
    
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
  
  print(myplot); 
  dev.off();
}
