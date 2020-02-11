# perform differential analysis
# default: all pair-wise comparison (A-B) + (B-C) + (A-C)
# custom: only compare two groups (A-C)
# time: compare consecutive groups (B-A) + (C-B)
# reference: all others against common reference (A-C) + (B-C)
# nested: (A-B)+(C-D) 
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param anal.type PARAM_DESCRIPTION, Default: 'default'
#' @param par1 PARAM_DESCRIPTION, Default: NULL
#' @param par2 PARAM_DESCRIPTION, Default: NULL
#' @param nested.opt PARAM_DESCRIPTION, Default: 'intonly'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformDEAnal
#' @export 
PerformDEAnal<-function (anal.type = "default", par1 = NULL, par2 = NULL, nested.opt = "intonly"){
  set.seed(1337);
  par1 <<- par1
  par2 <<- par2
  nested.opt <<-nested.opt
  myargs <- list()
  cls <- dataSet$cls
  dataSet$comp.type = anal.type
  grp.nms <- levels(cls)
  
  if (anal.type == "default") {
    inx = 0
    for (m in 1:(length(grp.nms) - 1)) {
      for (n in (m + 1):length(grp.nms)) {
        inx <- inx + 1
        myargs[[inx]] <- paste(grp.nms[m], "-", grp.nms[n], sep = "");
      }
    }
    filename = "SigGene_pairwise";
  } else if (anal.type == "time") {
    for (i in 2:length(grp.nms)) {
      myargs[[i - 1]] <- paste(grp.nms[i], "-", grp.nms[i-1], sep = "")
    }
    filename = "SigGene_time_series"
  } else if (anal.type == "custom") {
    grp.nms <- strsplit(par1, " vs. ")[[1]]
    myargs[[1]] <- paste(grp.nms, collapse = "-")
    dataSet$grp.nms = grp.nms;
    filename = paste("SigGene_", paste(grp.nms, collapse = "_vs_"), sep = "")
  } else if (anal.type == "reference") {
    ref <- par1;
    cntr.cls <- grp.nms[grp.nms != ref]
    myargs <- as.list(paste(cntr.cls, "-", ref, sep = ""));
    filename = paste("SigGene_reference_", ref, sep = "");
  } else if (anal.type == "nested") {
    grp.nms1 <- strsplit(par1, " vs. ")[[1]]
    grp.nms2 <- strsplit(par2, " vs. ")[[1]]
    if (all(grp.nms1 == grp.nms2)) {
      current.msg <<- paste("The two nested groups are the same. Please choose two different groups.")
      return(0)
    }
    grp.nms <- unique(c(grp.nms1, grp.nms2))
    if (nested.opt == "intonly") {
      myargs[[1]] <- paste("(", paste(grp.nms1, collapse = "-"), ")-(", paste(grp.nms2, collapse = "-"), ")", sep = "")
    } else {
      myargs[[1]] <- paste(grp.nms1, collapse = "-")
      myargs[[2]] <- paste(grp.nms2, collapse = "-")
      myargs[[3]] <- paste("(", paste(grp.nms1, collapse = "-"), ")-(", paste(grp.nms2, collapse = "-"), ")", sep = "")
    }
    filename = paste("SigGene_nested_", paste(paste(grp.nms1, collapse = "_vs_"), "_", paste(grp.nms2, collapse = "_vs_"), sep = ""), sep = "")
  } else {
    print(paste("Not supported: ", anal.type))
  }
  
  library(limma)
  design <- dataSet$design
  myargs[["levels"]] <- design
  contrast.matrix <- do.call(makeContrasts, myargs)
  if (dataSet$de.method == "limma") {
    if (is.null(dataSet$block)) {
      fit = lmFit(dataSet$data.norm, design)
    } else {
      corfit <- duplicateCorrelation(dataSet$data.norm, design, block = dataSet$block)
      fit <- lmFit(dataSet$data.norm, design, block = dataSet$block, correlation = corfit$consensus)
    }
    
    if (!is.fullrank(design)) {
      current.msg <<- paste("This metadata combination is not full rank! Please use other combination.")
      return(0)
    }
    
    df.residual <- fit$df.residual
    if (all(df.residual == 0)) {
      current.msg <<- paste("There is not enough replicates in each group (no residual degrees of freedom)!")
      return(0);
    }
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    topFeatures <- topTable(fit2, number = Inf, adjust.method = "fdr");
    
  } else if (dataSet$de.method == "deseq2"){
    # only for small data set (< 60)
    if(length(cls) > 60){
      current.msg <<- "For large sample size (>60), use limma or edgeR."; 
      return(0);
      
    }else{ # use microservice
      # use microservice
      print("Peforming DESeq2 ....");
      
      library(RSclient);
      rsc <- RS.connect();
      RS.assign(rsc, "my.dir", getwd()); 
      RS.eval(rsc, setwd(my.dir));
      
      dat.out <- list(data=dataSet, contrast.matrix = contrast.matrix);
      RS.assign(rsc, "dat.in", dat.out); 
      my.fun <- function(){
        library(DESeq2);
        dataSet <- dat.in$data;
        contrast.matrix <- dat.in$contrast.matrix;
        
        if (is.null(dataSet$sec.cls)){
          colData <- data.frame(dataSet$fst.cls)
          colnames(colData) <- "condition"
          dds <- DESeqDataSetFromMatrix(countData=round(dataSet$data.anot), colData = colData, design = ~condition);
        } else {
          colData <- data.frame(dataSet$fst.cls, dataSet$sec.cls, dataSet$cls);
          colnames(colData) <- c("condition", "type", "condition_type");
          dds <- DESeqDataSetFromMatrix(countData=round(dataSet$data.anot), colData = colData, design = ~condition_type);
        }   
        
        dds <- DESeq(dds, betaPrior=TRUE)
        vec <- as.numeric(c(0, contrast.matrix[,1]));
        res <- results(dds, contrast = vec, independentFiltering = FALSE, cooksCutoff = Inf);
        topFeatures <- data.frame(res@listData);
        rownames(topFeatures) <- rownames(res);
        nms <- colnames(topFeatures);
        nms[which(nms == "padj")] <- "adj.P.Val";
        nms[which(nms == "pvalue")] <- "P.Value";
        nms[which(nms == "log2FoldChange")] <- "logFC";
        colnames(topFeatures) <- nms;
        topFeatures <- topFeatures[c(2,1,3,4,5,6)];
        # order the result based on raw p
        ord.inx <- order(topFeatures$P.Value);
        topFeatures <- topFeatures[ord.inx, ];
        return(topFeatures);
      }
      RS.assign(rsc, my.fun);
      topFeatures <-  RS.eval(rsc, my.fun());
      RS.close(rsc);
    }
  } else {
    library(edgeR)
    y <- DGEList(counts = dataSet$data.anot, group = dataSet$cls)
    y <- calcNormFactors(y)
    y <- estimateGLMCommonDisp(y, design, verbose = FALSE)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, contrast = contrast.matrix)
    topFeatures <- topTags(lrt, n = Inf)$table
    nms <- colnames(topFeatures)
    nms[which(nms == "FDR")] <- "adj.P.Val"
    colnames(topFeatures) <- nms
  }
  dataSet$filename <- filename;
  dataSet$resTable <- topFeatures;
  dataSet <<- dataSet;
  return(1)
}
