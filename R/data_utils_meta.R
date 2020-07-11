##################################################
## R script for NetworkAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# read individual data from user, data.type: array
ReadIndData <- function(dataName, data.type){

    current.msg <<- "";
    meta.upload <<- FALSE; # upload data to resume
    dataSet <- .readTabData(dataName);

    # now need to remove low quality samples and genes
    data <- dataSet$data;
    meta.info <- dataSet$meta.info;

    smpl.num <- ncol(data);
    gene.num <- nrow(data);

    # remove smpls/exp with over half missing value
    good.inx<-apply(is.na(data), 2, sum)/nrow(data)<0.6;
    smpl.msg <- "";
    if(sum(!good.inx)>0){
        smpl.msg <- paste(sum(!good.inx), "low quality samples(>60% missing) removed.");
        print(smpl.msg);

        data <- data[,good.inx];
        if(ncol(data)/smpl.num < 0.5){
            current.msg <<- paste(smpl.msg, "Low quality data rejected!");;
            return("F");
        }

        # update meta information
        meta.info <- meta.info[good.inx, , drop=F];
    }

    if(ncol(data) < 4){
        current.msg <<- paste(smpl.msg, "The sample # (", ncol(data), ") is too small.");
        return("F");
    }

    # genes with 75% NA will be removed
    gd.inx<-apply(is.na(data), 1, sum)/ncol(data)<0.75;
    feat.msg <- "";
    if(sum(!gd.inx) > 0){
        data <- data[gd.inx,];
        feat.msg <- paste(sum(!gd.inx), "low quality genes (>75% missing) removed");
        if(nrow(data)/gene.num < 0.25){
            current.msg <<- paste(feat.msg, "Low quality data rejected.");
            return("F");
        }
        print(feat.msg);
    }

    if(nrow(data) < 10){ 
        current.msg <<- paste(feat.msg, "The gene# (", nrow(data), ") is too small (<10).");
        return("F");
    }

    # make an copy, only for testing different normalization
    dataSet$data.raw <- data;
    dataSet$data <- data;
    dataSet$type <- data.type;
    dataSet$meta.info <- meta.info;
    dataName <- dataSet$name;
    res <- RegisterData(dataSet);
    if(res == 1){
        return(dataName);
    }else{
        current.msg <<- paste("Cannot add data: ", dataName, ". ", current.msg, sep="");
        return("F");
    }
}

# read in the data and perform
# gene ID mapping using built in libraries
# matchMin is minimal matched probe (%)
# return the total matched gene number
ProcessIndData<-function(dataName, featureType, matchMin=0.5){

    dataSet <- readRDS(dataName);
    dataSet$id.type <- featureType;
    if(data.org != 'NA' & featureType != 'NA'){
        feature.vec <- rownames(dataSet$data.raw);
        minLvl <- length(feature.vec)*matchMin;
        entrez.id <-doAnnotation(feature.vec, featureType);

        hit.inx <- which(!is.na(entrez.id));
        matched.len <- length(hit.inx);
        if(matched.len > minLvl){
            res <- "Success!";
            data.orig <- dataSet$data.raw[hit.inx,];
            matched.entrez <- entrez.id[hit.inx];
            rownames(data.orig) <- matched.entrez;
            
            # now, deal with duplicated entrez id
            # first, average duplicate rows
            ave.data <- apply(data.orig, 2, myave, matched.entrez); 
            # then removed duplicated entries
            dup.inx <- duplicated(matched.entrez);
            int.mat <- ave.data[!dup.inx,];
            # update
            dataSet$data.orig <- int.mat;
            dataSet$id.type <- "entrez";
        }
        current.msg <<- paste("ID Conversion: ", "Total [", length(entrez.id), 
            "] Matched [", matched.len, "] Unmatched [", sum(is.na(entrez.id)),"]", collapse="\n"); 
    }else{ # no conversion will be performed
        dataSet$data.orig <- dataSet$data.raw;
        matched.len <- 9; # dummies
        minLvl <- 1;
    }
    RegisterData(dataSet);
    if(matched.len < minLvl){
        return(0);
    }
    return(matched.len);
}

# overwrite ave, => na.rm=T
myave <- function (x, ...) {
  n <- length(list(...))
  if (n) {
    g <- interaction(...)
    split(x, g) <- lapply(split(x, g), mean, na.rm=T)
  }
  else x[] <- FUN(x, na.rm=T)
  x
}

# read the uploaded data into memory
# return the meta-data information (multiple groups)
ReadDataForMetaInfo<-function(dataName){
    if(dataSet$name != dataName){
        dataSet <- readRDS(dataName);
    }
    return(colnames(dataSet$meta.info));
}

# here should first try to load the original data
# the data in the memory could be changed
GetGroupNames <- function(dataName){
    if(dataSet$name != dataName){
        dataSet <- readRDS(dataName);
    }
    levels(dataSet$cls);
}

GetDataDims <- function(dataName){
    if(dataSet$name != dataName){
        dataSet <- readRDS(dataName);
    }
    dm <- dim(dataSet$data);
    naNum <- sum(is.na(dataSet$data));
    return(c(dm, naNum));
} 

# obtain sample names and their class labels
GetSampleInfo <- function(dataName, clsLbl){
    if(dataSet$name != dataName){
        dataSet <- readRDS(dataName);
    }
    grpInfo <- dataSet$meta.info[[clsLbl]];
    grpLbls <- paste(levels(grpInfo), collapse="\n");
    smplInfo <- paste(Sample = colnames(dataSet$data.orig), "\t", Class=grpInfo, collapse="\n");
    return(c(grpLbls, smplInfo));
}

# given dataSet Name, sample name, and class name, do update
# note, for multiple #class, this set which one to use in the subsequent steps
# last one wins
UpdateSampleInfo<-function(dataName, clsLbl){
    
    print("updating sample info .... ");

    if(!exists("class.vec")){
        print("Could not find class label list!");
        return(0);
    }

    if(!exists("smpl.vec")){
        print("Could not find sample name list!");
        return(0);
    }

    if(length(class.vec) < 2){
        current.msg <<- "Add least two groups required!";
        return(0);
    }

    if(sum(class.vec != 'NA') < 2){
        current.msg <<- "Cannot be less than 2 groups";
        return(0);
    }

    if(dataSet$name != dataName){
        dataSet <- readRDS(dataName);
    }
    org.lvl.len <- length(levels(dataSet$meta.info[[clsLbl]]));
    if(org.lvl.len < length(class.vec)){
        current.msg <<- "You can not add new groups";
        return(0);
    }else if(org.lvl.len > length(class.vec)){
        current.msg <<- "To exclude a group, replace it with NA.";
        return(0);
    }

    # first update the meta info
    cls <- dataSet$meta.info[[clsLbl]];
    levels(cls) <- class.vec;

    data <- dataSet$data.orig;
    meta.info <- dataSet$meta.info;
 
    if(any(levels(cls) == 'NA')){
        rt.inx <- cls != 'NA';
        data <- data[,rt.inx];

        # also update the whole meta-info
        meta.info <- meta.info[rt.inx,,drop=FALSE];
        cls <- cls[rt.inx];
    }

    # need to re-construct the class, so that the level order  
    # are always alphabetic
    meta.info[[clsLbl]] <- factor(as.character(cls));

    # note, sample names could be removed (together with cls) as the whole row
    hit.inx <- colnames(data)%in%smpl.vec;
    dataSet$data.orig <- data[,hit.inx];

    # make sure the factor levels also dropped
    for(i in 1:length(meta.info)){
        meta.info[[i]] <- factor(meta.info[[i]][hit.inx]);
    }

    dataSet$meta.info <- meta.info;
    dataSet$cls <-  dataSet$meta.info[[clsLbl]];
    RegisterData(dataSet);
    gc();
    return(1);
}


# note, here also update data type array/count
PerformIndNormalization <- function(dataName, norm.opt, auto.opt, dataType){

    if(dataSet$name != dataName){
        dataSet <- readRDS(dataName);
    }
    msg <- NULL;
    data <- dataSet$data.orig;
    data <- PerformDataNormalization(data, norm.opt);
    if(length(data)==1 && data == 0){
        return(0);
    }
    msg <- paste(norm.msg, msg);

    if(auto.opt==1){
        row.nms <- rownames(data);
        col.nms <- colnames(data);
        data<-apply(data, 2, AutoNorm);
        msg <- paste(msg, "Autoscaling performed.", collapse=" ");
        rownames(data) <- row.nms;
        colnames(data) <- col.nms;
    }

    dataSet$data <- data;
    dataSet$type <- dataType;
    RegisterData(dataSet);
    current.msg <<- msg;
    return(1);
}

# normalize to zero mean and unit variance
AutoNorm<-function(x){
	(x - mean(x))/sd(x, na.rm=T);
}

######################################
## methods for merged expression data
#######################################

GlobalCutOff = list(
    logFC = 0,
    BHth = 0.05
)

# read meta-dataset previously processed
ReadMergedExpressTable <- function(dataName){
    current.msg <<- "";
    meta.upload <<- TRUE;
    dataSet <- .readTabData(dataName);
    common.matrix <- dataSet$data;
    meta.nms <- tolower(names(dataSet$meta.info));

    cls.inx <- grep("condition", meta.nms);
    if(length(cls.inx) == 0){
        current.msg <<- "No condition label found (#CLASS.condition)";
        return("F");
    }else{
        cls.inx <- cls.inx[1];
        cls.lbl <- dataSet$meta.info[[cls.inx]];
    }

    data.inx <- grep("dataset", meta.nms);
    if(length(data.inx) == 0){
        current.msg <<- "No dataset label found (#CLASS.dataset)";
        return("F");
    }else{
        data.inx <- data.inx[1];
        data.lbl <- dataSet$meta.info[[data.inx]];
        data.nms <- unique(as.character(data.lbl));
        
        # now create the mdata.all object
        mdata.all <- vector(mode="list", length=length(data.nms));
        names(mdata.all) <- data.nms;
        mdata.all <<- lapply(mdata.all, function(x){x=1});
    }

    if(length(grep("entrez.hsa", meta.nms)) > 0){
        data.org <<- "hsa"
        id.type <<- "entrez";
        shared.nms <- rownames(common.matrix);
        symbols <- doEntrez2SymbolMapping(shared.nms);
        names(symbols) <- shared.nms;
    }else if(length(grep("entrez.mmu", meta.nms)) > 0){
        data.org <<- "mmu"
        id.type <<- "entrez";
        shared.nms <- rownames(common.matrix);
        symbols <- doEntrez2SymbolMapping(shared.nms);
        names(symbols) <- shared.nms;
    }else{
        symbols <- NULL;
        inmex.org <<- "NA"
        id.type <<- 'NA';
    }
    data.org <<- unlist(strsplit( meta.nms[2], "[.]"))[3]
    inmex.meta.orig <- list(data=common.matrix,
                       id.type = id.type,
                       gene.symbls = symbols,
                       cls.lbl=factor(cls.lbl),
                       data.lbl=data.lbl);
    saveRDS(inmex.meta.orig, "inmex.meta.orig.rds");
    if(length(levels(as.factor(data.lbl))) == 1){
        return(2);
    }else{
        return(1);
    }    
}

GetMetaMetaInfo <- function(dataName){
    if(dataSet$name != dataName){
        dataSet <- readRDS(dataName);
    }
    return(colnames(dataSet$meta.info));
}

SetSelectedMetaMetaInfo <- function(dataName, meta0, meta1, block1){
    if(meta0 == "NA"){
        return(0);
    }else{
        if(dataSet$name != dataName){
            dataSet <- readRDS(dataName);
        }
         cls <- dataSet$meta.info[, meta0];
         block <- NULL;
         if(meta1 != "NA"){
            if(block1){
                block <- dataSet$meta.info[, meta1]
            }else{ # two factor
                cls <- interaction(dataSet$meta.info[, c(meta0, meta1)], sep = ".", lex.order = TRUE);
            }
         }
         dataSet$cls <- cls; # record main cls;
         dataSet$block <- block;
         RegisterData(dataSet);
         gc();
         return(levels(cls));
    }
}


# for gene-level meta-analysis
# function to set up results combining individual data analysis
# as well as to prepare for GO analysis
# no return, as set global 

SetupMetaStats <- function(BHth){
    
    GlobalCutOff$BHth <<- BHth;
    #all common genes
    inmex.meta <- readRDS("inmex_meta.rds");
    gene.ids <- rownames(inmex.meta$data);
    # meta.sig genes
    metade.genes <- rownames(meta.mat);

    # setup individual sig genes & stats
    # that overlap with meta.sig
    inmex.de <- list();

    pval.mat <- fc.mat <- matrix(nrow=nrow(meta.mat), ncol=sum(mdata.all==1));

    for(i in 1:length(inmex.ind)){
        de.res <- inmex.ind[[i]];

        hit.inx <- de.res[,2] <= BHth;
        hit.inx <- which(hit.inx); # need to get around NA
        inmex.de[[i]] <- rownames(de.res)[hit.inx];

        # only choose the genes that are also meta sig genes from in
        # individual analysis for display
        de.res <- de.res[metade.genes,];

        fc.mat[,i] <- de.res[,1];
        pval.mat[,i] <- de.res[,2];
    }
    names(inmex.de) <- names(inmex.ind);

    # calculate gain/loss
    deindst <- unique(unlist(inmex.de));
    gains=metade.genes[which(!(metade.genes %in% deindst))];
    losses=deindst[which(!(deindst %in% metade.genes))];
    all.de <- cbind(gene.ids %in% metade.genes, gene.ids %in% deindst);
    colnames(all.de) <- c("Meta-DE", "Individual-DE");
    vennC <- getVennCounts(all.de);
    if(inmex.meta$id.type == "entrez"){ 
        names(metade.genes) <- inmex.meta$gene.symbls[metade.genes];
        names(gains) <- inmex.meta$gene.symbls[gains];
        names(losses) <- inmex.meta$gene.symbls[losses];
    }

    # de genes from individual 
    de.len <- sapply(inmex.de, length);
    stat <- c(length(metade.genes), de.len);
    names(stat) <- c("Meta", substr(names(inmex.de), 0, nchar(names(inmex.de))-4));
    meta.stat <- list(
            stat = stat,
            de = metade.genes,
            idd = gains,
            loss = losses,
            venn = vennC
        );

    fc.mat <<- fc.mat;
    pval.mat <<- pval.mat;
    inmex.de <<- inmex.de;
    meta.stat <<- meta.stat;

    # save the result
    if(inmex.meta$id.type == "entrez"){ # row name gene symbols
        metade.nms <- inmex.meta$gene.symbls[metade.genes];
        res <- cbind(EntrezID=metade.genes, Name=metade.nms, meta.mat);
    }else{
        res <- cbind(ID=metade.genes, meta.mat);
    }
    write.csv(res, file=paste("meta_sig_genes_", inmex.method, ".csv", sep=""), row.names=F);
}

# perform differential analysis
# default: all pair-wise comparison (A-B) + (B-C) + (A-C)
# custom: only compare two groups (A-C)
# time: compare consecutive groups (B-A) + (C-B)
# reference: all others against common reference (A-C) + (B-C)
# nested: (A-B)+(C-D) 
PerformIndLimmaDE<-function(dataName, anal.type, par1=NULL, par2=NULL, nested.opt = "intonly"){
    if(dataSet$name != dataName){
        dataSet <- readRDS(dataName);
    }
    require(limma);
    current.msg <<- "";
    cls <- dSet$cls; 

    cls.lvls <- levels(cls);
    data <- dSet$data.norm;
    meta.info <- dSet$meta.info;

    design <- model.matrix(~ -1 + cls) # no intercept
    colnames(design) <- cls.lvls;

    if(is.null(dSet$block)){
        fit = lmFit(data, design);
    }else{
        # limma user guide chapter 8 P49 
        corfit<-duplicateCorrelation(data,design,block=dSet$block);
        fit = lmFit(data, design, block=dSet$block, correlation=corfit$consensus);
    }
    # sanity check
    if(!is.fullrank(design)){
         current.msg <<- paste("This metadata combination is not full rank! Please use other combination."); 
         return(0);
    }
    df.residual <- fit$df.residual;
    if (all(df.residual == 0)){
         current.msg <<- paste("There is not enough replicates in each group (no residual degrees of freedom)!"); 
         return(0);
    }

    myargs <- list();
    grp.nms <- cls.lvls;

    if(anal.type == 'default'){ # all pair-wise
        inx = 0;
        for(m in 1:(length(grp.nms)-1)){
            for(n in (m+1):length(grp.nms)){
                inx <- inx + 1;
                myargs[[inx]] <- paste(grp.nms[m], "-", grp.nms[n], sep="")
            }
        }
    }else if(anal.type == 'time'){
        for(i in 2:length(grp.nms)){
            myargs[[i-1]] <- paste(grp.nms[i], "-", grp.nms[i-1], sep="")
        }
    }else if(anal.type == 'custom'){
        grp.nms <- strsplit(par1, " vs. ")[[1]];
        myargs[[1]] <- paste(grp.nms, collapse="-");
    }else if(anal.type == 'reference'){
        ref <- par1;
        cntr.cls <- grp.nms[grp.nms!= ref];
        myargs <- as.list(paste(cntr.cls, "-", ref, sep=""));
    }else if(anal.type == 'nested'){
        grp.nms1 <- strsplit(par1, " vs. ")[[1]];
        grp.nms2 <- strsplit(par2, " vs. ")[[1]];
        if(all(grp.nms1 == grp.nms2)){
             current.msg <<- paste("The two nested groups are the same. Please choose two different groups."); 
             return(0);
        }
        grp.nms <- unique(c(grp.nms1, grp.nms2));
        if(nested.opt == "intonly"){
            myargs[[1]] <- paste("(", paste(grp.nms1, collapse="-"), ")-(", paste(grp.nms2, collapse="-"), ")", sep=""); 
        }else{
            myargs[[1]] <- paste(grp.nms1, collapse="-");
            myargs[[2]] <- paste(grp.nms2, collapse="-"); 
            myargs[[3]] <- paste("(", paste(grp.nms1, collapse="-"), ")-(", paste(grp.nms2, collapse="-"), ")", sep=""); 
        }
     }else{ # 
        print(paste('Not supported: ', anal.type));
     }

    myargs[["levels"]] <- design;
    contrast.matrix <- do.call(makeContrasts, myargs);
    fit2 <- contrasts.fit(fit, contrast.matrix);
    fit2 <- eBayes(fit2);

    resTable <- topTable(fit2, number=Inf, adjust.method="fdr");
    dSet$sig.orig <- resTable; # record orignal data for update
    anot.id <- rownames(resTable);
    gene.anot <- doEntrezIDAnot(anot.id);
    current.msg <<- current.msg; 
        
     # may need to update data and class
     hit.inx <- cls %in% grp.nms;
     if(sum(hit.inx) < length(hit.inx)){
         current.msg <<- paste(current.msg, "Only groups selected for comparisons: ", paste(grp.nms, collapse=", "), "are included.");
         cls <- factor(cls[hit.inx]);
         cls.lvls <- levels(cls);
         data <- data[,hit.inx];
         meta.info <- dSet$meta.info[hit.inx,];
     }

    dSet$data.stat <- data;
    dSet$cls.stat <- cls;
    dSet$meta.stat <- meta.info;
    dSet$sig.genes.anot <- gene.anot;
    RegisterData(dataSet);
    return (1);
}

GetProjectOrg <- function(){
    data.org;
}

GetProjectIDType <- function(){
    id.type;
} 

#  VENN DIAGRAM COUNTS AND PLOTS
getVennCounts <- function(x,include="both") {
	x <- as.matrix(x)
	include <- match.arg(include,c("both","up","down"))
	x <- sign(switch(include,
		both = abs(x),
		up = x > 0,
		down = x < 0
	))
	nprobes <- nrow(x)
	ncontrasts <- ncol(x)
	names <- colnames(x)
	if(is.null(names)) names <- paste("Group",1:ncontrasts)
	noutcomes <- 2^ncontrasts
	outcomes <- matrix(0,noutcomes,ncontrasts)
	colnames(outcomes) <- names
	for (j in 1:ncontrasts)
		outcomes[,j] <- rep(0:1,times=2^(j-1),each=2^(ncontrasts-j))
	xlist <- list()
	for (i in 1:ncontrasts) xlist[[i]] <- factor(x[,ncontrasts-i+1],levels=c(0,1))
	counts <- as.vector(table(xlist))
	structure(cbind(outcomes,Counts=counts),class="VennCounts")
}

#compute Cochran’s Q to help FEM/REM
# plot Q-Q plot for estimation
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

PlotDataProfile<-function(dataName, boxplotName, pcaName){
    dataSet <- readRDS(dataName);
    qc.boxplot2(dataSet$data, boxplotName);
    qc.pcaplot2(dataSet$data, pcaName);
}


qc.boxplot2 <- function(dat, imgNm){
    require('lattice');
    imgNm = paste(imgNm, "dpi", "72", ".png", sep="");
    subgene=10000;
    if (nrow(dat)>subgene) {
        set.seed(28051968);
        sg  = sample(nrow(dat), subgene)
        Mss = dat[sg,,drop=FALSE]
    } else {
        Mss = dat
    }

    subsmpl=100;
    if (ncol(Mss)>subsmpl) {
        set.seed(28051968);
        ss  = sample(ncol(Mss), subsmpl)
        Mss = Mss[,ss,drop=FALSE]
    } else {
        Mss = Mss
    }

    sample_id = rep(seq_len(ncol(Mss)), each = nrow(Mss));
    values  = as.numeric(Mss)
    formula = sample_id ~ values

    box = bwplot(formula, groups = sample_id, layout = c(1,1), as.table = TRUE,
        strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
        horizontal = TRUE,
        pch = "|",  col = "black", do.out = FALSE, box.ratio = 2,
        xlab = "", ylab = "Samples",
        fill = "#1c61b6AA",
        panel = panel.superpose,
        scales = list(x=list(relation="free"), y=list(axs="i")),
        ylim = c(ncol(Mss)+0.7,0.3),
        prepanel = function(x, y) {
          list(xlim = quantile(x, probs = c(0.01, 0.99), na.rm=TRUE))
        },
        panel.groups = function(x, y, ...) {
          panel.bwplot(x, y, ...)
        })

    Cairo(file=imgNm, width=460, height=420, type="png", bg="white");
    print(box);
    dev.off();
}

qc.pcaplot2 <- function(x, imgNm){
    imgNm = paste(imgNm, "dpi", "72", ".png", sep="");
    require('lattice');
    pca <- prcomp(t(na.omit(x)));
    names <- colnames(x);
    pca.res <- as.data.frame(pca$x);
    # increase xlim ylim for text label
    xlim <- GetExtendRange(pca.res$PC1);
    ylim <- GetExtendRange(pca.res$PC2);
    pcafig = xyplot(PC2~PC1, data=pca.res, pch=19, cex=1,aspect = "iso", xlim = xlim, ylim=ylim,
        panel=function(x, y, ...) {
               panel.xyplot(x, y, ...);
               ltext(x=x, y=y, labels=names, pos=1, offset=1, cex=0.8, col="magenta")
            })

    Cairo(file=imgNm, width=480, height=480, type="png", bg="white");
    print(pcafig);
    dev.off();
}

PlotMetaPCA <- function(imgNm, dpi, format,factor){
  inmex.meta <- readRDS("inmex_meta.rds");
  x=inmex.meta[["data"]]
  dpi = as.numeric(dpi);
  imgNm = paste(imgNm, "dpi", dpi, ".", format, sep="");
  library('lattice');
  library('ggplot2');
  pca <- prcomp(t(na.omit(x)));
  imp.pca<-summary(pca)$importance;
  xlabel = paste0("PC1"," (", 100*round(imp.pca[2,][1], 3), "%)")
  ylabel = paste0("PC2"," (", 100*round(imp.pca[2,][2], 3), "%)")
  names <- colnames(x);
  pca.res <- as.data.frame(pca$x);
  # increase xlim ylim for text label
  xlim <- GetExtendRange(pca.res$PC1);
  ylim <- GetExtendRange(pca.res$PC2);
  if(factor != "NA"){
    #Factor = as.vector(dataSet$meta.info[,factor])
  }else{
    #Factor = dataSet$meta.info[,1];
  }
  Conditions = factor(inmex.meta$cls.lbl)
  Datasets = factor(inmex.meta$data.lbl)
  pcafig = ggplot(pca.res, aes(x=PC1, y=PC2,  color=Conditions ,shape=Datasets)) +
    geom_point(size=4, alpha=0.5) + xlim(xlim)+ ylim(ylim) + xlab(xlabel) + ylab(ylabel);
  
  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", unit="in", dpi=dpi);
  print(pcafig);
  dev.off();
  
}


GetMetaSummary<- function(){
   inmex.meta <- readRDS("inmex_meta.rds");
    sel.nms <- unique(inmex.meta$data.lbl)
    sel.nms <- paste(sel.nms, collapse="; ")
    cls.lbls <- unique(inmex.meta$cls.lbl)
    cls.lbls <- paste(cls.lbls, collapse="; ")
    return(c(length(colnames(inmex.meta$data)),nrow(inmex.meta$data), sel.nms, cls.lbls))
}

GetMetaDatasets<- function(){
  sel.nms <- names(mdata.all)[mdata.all==1];
  return(sel.nms);
}

SetSelMetaData<- function(selNm){
    selDataNm <<- selNm
}

# retrun the json obj
SaveMetaClusterJSON <- function(fileName, clustOpt, opt){

    inmex.meta <- readRDS("inmex_meta.rds");
    datanm.vec <- names(mdata.all)[mdata.all==1];
    dat.inx <- inmex.meta$data.lbl %in% datanm.vec;
    dat <- inmex.meta$data[, dat.inx, drop=F]; 

    # need to deal with missing values 
    dat <- na.omit(dat);

    pca3d <- list();
    if(clustOpt == "pca"){
        if(opt == "all"){
            pca <- prcomp(t(dat), center=T, scale=T);
            }else{
            dat = dat[which(rownames(dat) %in% loadEntrez),]
            pca <- prcomp(t(dat), center=T, scale=T);
            }

        imp.pca<-summary(pca)$importance;
        pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
        coords <- data.frame(t(signif(pca$x[,1:3], 5)));
    }else{
        require('Rtsne');
        ndat <- as.matrix(t(dat));
        max.perx <- floor((nrow(ndat)-1)/3);
        if(max.perx > 30){
            max.perx <- 30;
        }
        res <- Rtsne(ndat, dims = 3, perplexity=max.perx);
        pca3d$score$axis <- paste("t-SNE dim ", 1:3, sep="");
        coords <- data.frame(t(signif(res$Y, 5)));
    }

    colnames(coords) <- NULL; 
    pca3d$score$xyz <- coords;
    pca3d$score$name <- colnames(dat);

    facA <- as.character(inmex.meta$cls.lbl[dat.inx]);
    if(all.numeric(facA)){
        facA <- paste("Group", facA);
    }
    pca3d$score$facA <- facA;

    facB <-  as.character(inmex.meta$data.lbl[dat.inx]);
    if(all.numeric(facB)){
        facB <- paste("Group", facB);
    }
    pca3d$score$facB <- facB;

    # now set color for each group
    cols <- unique(GetColorSchema(facB));
    rgbcols <- col2rgb(cols);
    cols <- apply(rgbcols, 2, function(x){paste("rgba(", paste(x, collapse=","), ",1)", sep="")});
    pca3d$score$colors <- cols;

    # add shape sphere, triangles, square, pentagon (first two)
    pca3d$score$shapes <- c("sphere", "triangle");

    mypos <- t(coords);
    colnames(mypos) <- paste("Dim", 1:3, sep="");
    coords <- data.frame(Class=facA, Data=facB, mypos);
    write.csv(coords, file="networkanalyst_3d_pos.csv");

    require(RJSONIO);
    pca3d$org = data.org
    pca3d$analType=anal.type
    pca3d$naviString = "PCA 3D"
    dataSet$jsonNms$pcascore <<- fileName
    json.mat <- toJSON(pca3d, .na='null');
    partialToBeSaved <<- c(partialToBeSaved, c(fileName))
    sink(fileName);
    cat(json.mat);
    sink();
    current.msg <<- "Annotated data is now ready for 3D visualization!";
    return(1);
}


SaveMetaClusterLoadingJSON <- function(fileName, clustOpt, nb){

    inmex.meta <- readRDS("inmex_meta.rds");
    datanm.vec <- names(mdata.all)[mdata.all==1];
    nb = as.numeric(nb)
    dat.inx <- inmex.meta$data.lbl %in% datanm.vec;
    dat <- inmex.meta$data[, dat.inx, drop=F]; 

    # need to deal with missing values 
    dat <- na.omit(dat);
    variances = apply(dat,1, function(x){var(x)})
    df = data.frame(var = variances, inx = seq.int(1,length(variances)))
    df = df[order(-df$var),]
    inx = df$inx[c(1:nb)]
    dat = dat[inx,];

    pca3d <- list();

    pca <- prcomp(t(dat), center=T, scale=T);    
        imp.pca<-summary(pca)$importance;
        pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
        coords <- data.frame(t(signif(pca$rotation[,1:3], 5)));
        
    colnames(coords) <- NULL; 
    pca3d$score$xyz <- coords;
    pca3d$score$name <- doEntrez2SymbolMapping(rownames(pca$rotation));
    pca3d$score$entrez <- rownames(pca$rotation);

    loadEntrez <<- pca3d$score$entrez
    mypos <- t(coords);
    colnames(mypos) <- paste("Dim", 1:3, sep="");
    
    coords <- data.frame(mypos);
    write.csv(coords, file="networkanalyst_loadings_3d_pos.csv");
   
    require(RJSONIO);
    partialToBeSaved <<- c(partialToBeSaved, c(fileName))
    dataSet$jsonNms$pcaload <<- fileName
    json.mat <- toJSON(pca3d, .na='null');
    sink(fileName);
    cat(json.mat);
    sink();
    current.msg <<- "Annotated data is now ready for 3D visualization!";
    return(1);
}

GetDatasetNamesString <- function(){
    inmex.meta <- readRDS("inmex_meta.rds");
    paste(unique(inmex.meta$data.lbl), collapse="||");
}

PerformBatchCorrection <- function(){
    .prepare.batch();
    .perform.computing();
    # no need to save, already done
}

.prepare.batch<-function(){
    my.fun <- function(){
        library('sva');
        inmex.meta <- readRDS("inmex_meta.rds");
        data.lbl <- inmex.meta$data.lbl
        pheno <- data.frame(cbind(inmex.meta$cls.lbl, data.lbl));
        modcombat <- model.matrix(~1, data=pheno)
        batch <- data.lbl;
        inmex.meta$data <- ComBat(dat=inmex.meta$data, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE);
        saveRDS(inmex.meta, "inmex_meta.rds");
    }
    dat.in <- list(my.fun=my.fun);
    saveRDS(dat.in, file="dat.in.rds");
    return(1)
}
