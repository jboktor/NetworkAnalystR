# determine if all annotated data are ready for meta-analysis
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[stringi]{stri_isempty}},\code{\link[stringi]{stri_length}},\code{\link[stringi]{stri_numbytes}},\code{\link[stringi]{stri_width}}
#'  \code{\link[ipred]{print.classbagg}}
#'  \code{\link[globals]{Globals}}
#'  \code{\link[gdata]{mapLevels}}
#'  \code{\link[BiocGenerics]{row+colnames}},\code{\link[BiocGenerics]{paste}},\code{\link[BiocGenerics]{sets}},\code{\link[BiocGenerics]{unlist}},\code{\link[BiocGenerics]{lapply}},\code{\link[BiocGenerics]{t}},\code{\link[BiocGenerics]{cbind}},\code{\link[BiocGenerics]{unique}},\code{\link[BiocGenerics]{order}}
#'  \code{\link[AnnotationDbi]{Bimap-toTable}}
#'  \code{\link[S4Vectors]{Vector-class}}
#'  \code{\link[h2o]{scale}}
#'  \code{\link[matrixStats]{matrixStats-methods}}
#'  \code{\link[dplyr]{reexports}}
#'  \code{\link[utils]{write.table}}
#'  \code{\link[fields]{compactToMat}}
#' @rdname CheckMetaDataIntegrity
#' @export 
#' @importFrom stringi length
#' @importFrom ipred print
#' @importFrom globals names
#' @importFrom gdata factor levels
#' @importFrom BiocGenerics rownames paste intersect unlist lapply t cbind unique order
#' @importFrom AnnotationDbi colnames ncol nrow
#' @importFrom S4Vectors rep c
#' @importFrom h2o scale
#' @importFrom matrixStats vector
#' @importFrom dplyr data.frame
#' @importFrom utils write.table
#' @importFrom fields matrix
CheckMetaDataIntegrity<-function(){
  
  performedDE <<- FALSE;
  
  if(stringi::length(mdata.all) == 0){
    current.msg <<-"Please upload your data or try our example datasets!";
    ipred::print(current.msg);
    return(0);
  }
  
  include.inx <- mdata.all==1;
  if(sum(include.inx) == 0){
    current.msg <<-"No dataset is selected for analysis!";
    ipred::print(current.msg);
    return(0);
  }
  
  sel.nms <- globals::names(mdata.all)[include.inx];
  clss <- list();
  if(meta.upload){
    # update meta data only for select/deselect datasets
    inmex.meta.orig <- readRDS("inmex.meta.orig.rds");
    hit.inx <- inmex.meta.orig$data.lbl %in% sel.nms;
    data <- inmex.meta.orig$data[, hit.inx];
    id.type <- inmex.meta.orig$id.type;
    cls.lbl <- gdata::factor(inmex.meta.orig$cls.lbl[hit.inx]);
    data.lbl <- gdata::factor(inmex.meta.orig$data.lbl[hit.inx]);
    common.matrix <- data;
    nms.vec <<- BiocGenerics::rownames(inmex.meta.orig$data);
    smps.vec <<- data.lbl
  }else{   
    # first check that all class labels are consistent
    dataSet <- readRDS(sel.nms[1]);
    lvls <- gdata::levels(dataSet$cls);
    id.type <- dataSet$id.type;
    clss[[1]] = dataSet$cls;
    nms <- BiocGenerics::rownames(dataSet$data);
    shared.nms <- nms;
    for(i in 2:stringi::length(sel.nms)){
      dataSet <- readRDS(sel.nms[i]);
      clss[[i]] = dataSet$cls;
      # check if class label is consistent
      if(!all(gdata::levels(dataSet$cls) == lvls)){
        current.msg <<- BiocGenerics::paste(sel.nms[i], "has different group labels", 
                              BiocGenerics::paste(gdata::levels(dataSet$cls), collapse=":"), 
                              "from", sel.nms[1], BiocGenerics::paste(lvls, collapse=":"));
        ipred::print(current.msg);
        return(0);
      }
      
      # check and record if there is common genes            
      shared.nms <- BiocGenerics::intersect(shared.nms, BiocGenerics::rownames(dataSet$data));
      if(stringi::length(shared.nms) < 10){
        current.msg <<- BiocGenerics::paste(sel.nms[i], "has less than 10 common genes/probes from previous data sets");
        ipred::print(current.msg);
        return(0);
      }
      
      # check gene id type
      if(dataSet$id.type != id.type){
        current.msg <<- BiocGenerics::paste(sel.nms[i], "has different gene/probe ID from", sel.nms[1]);
        ipred::print(current.msg);
        return(0);
      }
    }
    
    nrepu<<-BiocGenerics::unlist(BiocGenerics::lapply(clss,FUN=function(x)stringi::length(x)));      
    
    ipred::print("Passed exp condition check!");
    
    # now construct a common matrix to faciliated plotting across all studies
    dataName <- sel.nms[1];
    dataSet <- readRDS(dataName);
    common.matrix <- dataSet$data[shared.nms, ];
    nms.vec = BiocGenerics::rownames(dataSet$data);
    smps.vec = AnnotationDbi::colnames(dataSet$data);
    data.lbl <- S4Vectors::rep(dataName, AnnotationDbi::ncol(common.matrix));
    cls.lbl <- dataSet$cls;
    
    for(i in 2:stringi::length(sel.nms)){
      dataName <- sel.nms[i];
      dataSet <- readRDS(dataName);
      ndat <- dataSet$data[shared.nms, ];
      nms.vec = S4Vectors::c(nms.vec, BiocGenerics::rownames(dataSet$data));
      smps.vec = S4Vectors::c(smps.vec, AnnotationDbi::colnames(dataSet$data));
      plot.ndat <- BiocGenerics::t(h2o::scale(BiocGenerics::t(ndat)));
      common.matrix <- BiocGenerics::cbind(common.matrix, ndat);
      data.lbl <- S4Vectors::c(data.lbl, S4Vectors::rep(dataName, AnnotationDbi::ncol(dataSet$data[,])));
      cls.lbl <- S4Vectors::c(cls.lbl, dataSet$cls);
    }
    cls.lbl <- gdata::factor(cls.lbl);
    gdata::levels(cls.lbl) <- lvls;
    
    smps.nms = AnnotationDbi::colnames(common.matrix)
    if(stringi::length(BiocGenerics::unique(smps.nms)) != stringi::length(smps.nms)){
      data.nb = stringi::length(BiocGenerics::unique(data.lbl));
      data.vec = matrixStats::vector()
      for(i in 1:data.nb){
        data.vec[i] = paste0("d", i);
      }
      gdata::levels(data.lbl) = data.vec;
      AnnotationDbi::colnames(common.matrix) = make.unique(BiocGenerics::paste(data.vec, smps.nms, sep="_"));
      
      dataSet <- readRDS(sel.nms[1]);
      AnnotationDbi::colnames(dataSet$data) = BiocGenerics::paste("d1", AnnotationDbi::colnames(dataSet$data), sep="_");
      RegisterData(dataSet);
      
      for(i in 2:stringi::length(sel.nms)){
        dataSet <- readRDS(sel.nms[i]);
        AnnotationDbi::colnames(dataSet$data) = paste0("d",i,"_",AnnotationDbi::colnames(dataSet$data));
        # check if class label is consistent
        RegisterData(dataSet);
      }
      smps.vec <<- smps.nms;
      current.msg <<- BiocGenerics::paste("Duplicated sample names detected, samples have been renamed to make them unique.");
    }
    
    # note: index by entrez, gene symbol DO have duplicates
    BiocGenerics::rownames(common.matrix) <- shared.nms;
    
    # resort data, first on data.lbl, then on class lbl
    ord.inx <- BiocGenerics::order(data.lbl, cls.lbl);
    common.matrix <- data.matrix(common.matrix[,ord.inx]);
    cls.lbl <- cls.lbl[ord.inx];
    data.lbl <- data.lbl[ord.inx];
    smps.vec <- smps.vec[ord.inx];
    nms.vec = BiocGenerics::unique(nms.vec)
    nms.vec <<- nms.vec
    smps.vec <<- smps.vec
  }
  
  if(AnnotationDbi::ncol(common.matrix) > 1000){  # max sample number allow 1000
    current.msg <<- BiocGenerics::paste("Total combined sample #:", AnnotationDbi::ncol(common.matrix), "(exceed the limit: 1000!)");
    return(0);
  }
  
  # save the meta-dataset
  res <- dplyr::data.frame(AnnotationDbi::colnames(common.matrix), cls.lbl, data.lbl, BiocGenerics::t(common.matrix));
  AnnotationDbi::colnames(res) <- S4Vectors::c('#NAME', '#CLASS.condition', BiocGenerics::paste('#CLASS.dataset',id.type, data.org, sep="."), BiocGenerics::rownames(common.matrix));
  utils::write.table(BiocGenerics::t(res), sep="\t", file="NetworkAnalyst_merged_data.txt", col.names=F, quote=FALSE);
  
  # need to set up the data for plotting (boxplot, heatmpa) so 
  # we need to scale row for each dataset in order to elimiate the maganitude difference 
  plot.matrix <- fields::matrix(NA, nrow=AnnotationDbi::nrow(common.matrix), ncol=AnnotationDbi::ncol(common.matrix));
  BiocGenerics::rownames(plot.matrix) <- BiocGenerics::rownames(common.matrix);
  AnnotationDbi::colnames(plot.matrix) <- AnnotationDbi::colnames(common.matrix);
  for(i in 1:stringi::length(sel.nms)){
    data.inx <- data.lbl == sel.nms[i];
    plot.matrix[,data.inx] <- BiocGenerics::t(h2o::scale(BiocGenerics::t(common.matrix[,data.inx])));
  }
  
  # if entrez, get symbols for display
  shared.nms <- BiocGenerics::rownames(common.matrix);
  if(id.type == "entrez"){ 
    symbols <- doEntrez2SymbolMapping(shared.nms);
  }else{ # display itself
    symbols <- shared.nms;
  }
  globals::names(symbols) <- shared.nms;
  
  inmex.meta <- list(data=common.matrix,
                     plot.data=plot.matrix,
                     id.type = id.type,
                     gene.symbls = symbols,
                     cls.lbl=gdata::factor(cls.lbl),
                     data.lbl=data.lbl);
  
  saveRDS(inmex.meta, "inmex_meta.rds");
  smps.vec <<- AnnotationDbi::colnames(common.matrix);
  
  # setup common stats gene number, smpl number, grp info
  current.msg <<- BiocGenerics::paste("Sample #:", AnnotationDbi::ncol(inmex.meta$data),
                        "Common ID #:", AnnotationDbi::nrow(inmex.meta$data), 
                        "Condition:", BiocGenerics::paste(gdata::levels(inmex.meta$cls.lbl), collapse=" vs. "));
  return(1);
}
