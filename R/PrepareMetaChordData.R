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
#' @rdname PrepareMetaChordData
#' @export 
PrepareMetaChordData <- function(){
  BHth <- GlobalCutOff$BHth;
  logFC <- GlobalCutOff$logFC;
  sel.nms <- names(mdata.all)[mdata.all==1];
  row.nm <- length(sel.nms);
  
  # update inmex.ind for only selected dataset
  hit.inx <-  names(inmex.ind) %in% sel.nms;
  inmex.ind <- inmex.ind[hit.inx];
  dataNms <- names(inmex.ind);
  newNms <- substring(dataNms,0, nchar(dataNms)-4);
  
  if(meta.selected){
    newNms <- c(newNms, "meta");
    row.nm <- row.nm + 1;
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  hit.mat <- matrix(0, nrow=row.nm, ncol=nrow(inmex.meta$data));
  nm.mat <- matrix(NA, nrow=row.nm, ncol=nrow(inmex.meta$data));
  shared.ids <- rownames(inmex.meta$data);
  shared.sbls <- inmex.meta$gene.symbls;
  
  colnames(hit.mat) <- colnames(nm.mat) <- shared.sbls;
  rownames(hit.mat) <- rownames(nm.mat) <- newNms;
  
  chord.vec.nms <- newNms;
  chord.res <- list();
  tot.count <- 0; # need to see how many chords
  for(i in 1:length(inmex.ind)){
    dataSet = readRDS(sel.nms[i]);
    nm <- newNms[i];
    res.mat <- inmex.ind[[i]][shared.ids, ];
    hit.inx <- abs(res.mat[,1]) >= logFC & res.mat[,2] <= dataSet$pval;
    hit.mat[i, hit.inx] <- 1;
    nm.mat[i, hit.inx] <- nm;
    chord.res[[nm]] <- shared.sbls[hit.inx];
    tot.count <- tot.count + sum(hit.inx);
  }
  
  if(meta.selected){# add meta
    i  <- length(inmex.ind) + 1;
    hit.inx <- shared.ids %in% meta.stat$de;
    hit.mat[i, hit.inx] <- 1;
    nm.mat[i, hit.inx] <- newNms[i];
    chord.res[[newNms[i]]] <- shared.sbls[hit.inx];
    tot.count <- tot.count + sum(hit.inx);
  }
  
  if(tot.count > 2000){
    current.msg <<- paste("Chord diagrams is effective to display relationships for less than 1000 genes. The results contain", tot.count, 
                          "of genes (max. allowed: 2000). You can try Venn diagram or adjust threshold to select most significant genes.")
    print(current.msg);
    return(NULL);
  }
  
  # keep gene with at least 1 hit across all datasets
  keep.inx <- apply(hit.mat, 2, sum) > 0;
  hit.mat <- hit.mat[,keep.inx];
  nm.mat <- nm.mat[,keep.inx];
  shared.sbls <- shared.sbls[keep.inx];
  
  ##############Links#######
  #links & arcs data 
  links <- arcs <- list();
  
  # arcs data
  de.num <- apply(hit.mat, 1, sum);
  de.prct <- as.numeric(de.num/sum(de.num));
  
  for(i in 1:length(newNms)){
    arcs[[i]] <- list(
      name = newNms[i],
      label = newNms[i],
      size = as.numeric(de.num[i]),
      val = de.prct[i]
    )
  }
  
  # now, get unique names by each inserting in the gene names
  lnk.genes <- colnames(nm.mat);
  nm.mat <- apply(nm.mat, 1, 
                  function(x){
                    gd.hits <- !is.na(x); 
                    x[gd.hits] <-paste(x[gd.hits], "*", lnk.genes[gd.hits], sep="");
                    x;
                  });
  
  nm.mat <- unname(t(nm.mat));
  for(l in 1:nrow(nm.mat)){
    x <- nm.mat[l,];      
    hit.inx <- !is.na(x);
    links$name <- c(links$name, x[hit.inx]);
    
    y.mat <- nm.mat[,hit.inx,drop=F];
    new.mat <- apply(y.mat, 2, 
                     function(d){
                       myd <- d[-l];
                       if(all(is.na(myd))){ # link to itself?
                         d[l];
                       }else{
                         myd[!is.na(myd)];
                       }
                     });
    links$imports <- c(links$imports,new.mat);
  }
  
  # need to reformat
  links.new <- list();
  for(l in 1:length(links$name)){
    impt<-links$imports[[l]];
    if(length(impt) == 1){
      impt <- as.matrix(impt);
    }
    links.new[[l]] <- list(
      name=links$name[l],
      imports=impt
    )
  }
  
  ###################
  #colors & labels
  require(RColorBrewer);  
  if(length(newNms) > 9){
    colors <- brewer.pal(length(newNms),"Set3");
  }else if (length(newNms) > 2){
    colors <- brewer.pal(length(newNms),"Dark2");
  }else {
    colors <- c("#7570B3", "#D95F02");
  }
  labels <- rep("", length(newNms));
  names(labels) <- names(colors) <- newNms;
  
  ###################
  # links done!, now get pvalues/data for every pairs of chord
  dataNms <- rownames(hit.mat);
  geneNms <- colnames(hit.mat);
  entrezIDs <- names(shared.sbls);
  weights <- data <- list();
  
  # also need to setup lookup json
  lookup <- list();
  for(i in 1:length(dataNms)){
    nlst <- vector(length = length(dataNms)-1, mode = "list"); 
    nm <- dataNms[i];
    names(nlst) <- dataNms[-i]
    lookup[[nm]] <- nlst;
  }
  
  item.count <- 0;
  for(i in 1:ncol(nm.mat)){
    geneNm <- geneNms[i]; 
    entrez <- entrezIDs[i];
    hit.inx <- !is.na(nm.mat[,i]);
    datNms <- dataNms[hit.inx];
    
    if(length(datNms) == 1){
      item.count <- item.count + 1;
      id <-paste(datNms,"*",datNms,"*",geneNm, sep="");
      weights[[id]]<-1;
      lookup[[datNms]][[datNms]][[geneNm]] <- list(
        entrez = entrez,
        genesymbol = geneNm
      )
      data[[item.count]]<-list(
        fromModule = datNms,
        toModule= datNms,
        arc=geneNm,
        fromColor=as.character(colors[datNms]),
        toColor=as.character(colors[datNms])
      );
    }else{
      for(m in 1:(length(datNms)-1)){
        nm1 <- datNms[m];
        for(n in (m+1):length(datNms)){
          item.count <- item.count + 1;
          nm2 <- datNms[n];
          id1<- paste(nm1,"*",nm2,"*",geneNm, sep="");
          id2<- paste(nm2,"*",nm1,"*",geneNm, sep="");
          weights[[id1]]<-weights[[id2]]<-1;
          
          lookup[[nm1]][[nm2]][[geneNm]] <- list(
            entrez = entrez,
            genesymbol = geneNm
          )
          lookup[[nm2]][[nm1]][[geneNm]] <- list(
            entrez = entrez,
            genesymbol = geneNm
          )
          # only from large module (src) to small (no need for reverse)
          if(de.num[nm1] >= de.num[nm2]){
            data[[item.count]]<-list(
              fromModule = nm1,
              toModule= nm2,
              arc=geneNm,
              fromColor=as.character(colors[nm1]),
              toColor=as.character(colors[nm2])
            );
          }else{
            data[[item.count]]<-list(
              fromModule = nm2,
              toModule=nm1,
              arc=geneNm,
              fromColor=as.character(colors[nm2]),
              toColor=as.character(colors[nm1])
            );
          }
        }
      }
    }
  }
  
  # need to remove 3rd level list name 
  for(m in 1:length(dataNms)){ # need to remove itself
    nm1 <- dataNms[m];
    for(n in m:length(dataNms)){
      nm2 <- dataNms[n];
      nlist <- lookup[[nm1]][[nm2]];
      lookup[[nm1]][[nm2]] <- lookup[[nm2]][[nm1]] <- unname(nlist);
    }
  }
  
  chordData <- list(
    arcs = arcs,
    links = links.new,
    data = data,
    weights = weights,
    colors = colors,
    labels = labels
  )
  
  chord.vec.nms <<- chord.vec.nms;
  chord.res <<- chord.res;
  
  return (
    list(
      chordData = chordData,
      lookup = lookup
    )
  );
}
