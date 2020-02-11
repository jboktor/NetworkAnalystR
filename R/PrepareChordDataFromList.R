#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param newDat PARAM_DESCRIPTION
#' @param uniq.enIDs PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PrepareChordDataFromList
#' @export 
PrepareChordDataFromList <- function(newDat, uniq.enIDs){
  # now combine lists into a single matrix
  all.nms <- listNms
  nm.used <- NULL;
  for(i in 1:length(all.nms)){
    dataNm <- all.nms[i];
    if(mdata.all[[dataNm]] == 1){ # selected
      nm.used <- c(nm.used, dataNm);
    }
  }
  if(anal.type == "metadata" & meta.selected){
    nm.used <- c(nm.used, "meta");
  }
  hit.mat <- fc.mat <- matrix(0, nrow=length(nm.used), ncol=length(uniq.enIDs));
  nm.mat <- matrix(NA, nrow=length(nm.used), ncol=length(uniq.enIDs));
  shared.ids <- uniq.enIDs;
  shared.sbls <- doEntrez2SymbolMapping(uniq.enIDs);
  
  dataNms <-nm.used;
  colnames(hit.mat) <- colnames(fc.mat) <- colnames(nm.mat) <- shared.sbls;
  rownames(hit.mat) <- rownames(fc.mat) <- rownames(nm.mat) <- dataNms;
  names(shared.sbls) <- shared.ids;
  
  # push to parent env.
  shared.sbls <<- shared.sbls;
  
  for(i in 1:length(dataNms)){
    nm <- dataNms[i];
    exp.vec <- newDat[[nm]];
    hit.inx <- shared.ids %in% names(exp.vec);
    hit.mat[i, hit.inx] <- 1;
    nm.mat[i, hit.inx] <- nm;
    fc.mat[i, hit.inx] <- as.numeric(exp.vec);
  }
  
  ##############Links#######
  ## links & arcs data 
  ###########################
  
  links <- arcs <- list();
  
  # arcs data
  de.num <- apply(hit.mat, 1, sum);
  de.prct <- as.numeric(de.num/sum(de.num));
  
  for(i in 1:length(dataNms)){
    arcs[[i]] <- list(
      name = dataNms[i],
      label = dataNms[i],
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
  if(length(dataNms) > 9){
    colors <- brewer.pal(length(dataNms),"Set3");
  }else if (length(dataNms) > 2){
    colors <- brewer.pal(length(dataNms),"Dark2");
  }else {
    colors <- c("#7570B3", "#D95F02");
  }
  labels <- rep("", length(dataNms));
  names(labels) <- names(colors) <- dataNms;
  
  ###################
  # links done!, now get weights/pvalues/data for every pairs of chord
  dataNms <- rownames(fc.mat);
  geneNms <- colnames(fc.mat);
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
    fc.vals <- as.numeric(fc.mat[hit.inx,i]);
    
    datNms <- dataNms[hit.inx];
    
    if(length(datNms) == 1){
      item.count <- item.count + 1;
      id <-paste(datNms,"*",datNms,"*",geneNm, sep="")
      #weights[[id]]<-fc.vals;
      weights[[id]]<-1;
      lookup[[datNms]][[datNms]][[geneNm]] <- list(
        entrez = entrez,
        genesymbol = geneNm,
        weight = fc.vals
      )
      data[[item.count]]<-list(
        fromModule = datNms,
        toModule= datNms,
        arc=geneNm,
        weight=fc.vals,
        fromColor=as.character(colors[datNms]),
        toColor=as.character(colors[datNms])
      );
    }else{
      for(m in 1:(length(datNms)-1)){
        nm1 <- datNms[m];
        fc1 <- fc.vals[m];
        for(n in (m+1):length(datNms)){
          item.count <- item.count + 1;
          nm2 <- datNms[n];
          fc2 <- fc.vals[n];
          id1<- paste(nm1,"*",nm2,"*",geneNm, sep="");
          id2<- paste(nm2,"*",nm1,"*",geneNm, sep="");
          wt <- 1;
          weights[[id1]]<-weights[[id2]]<-1;
          
          lookup[[nm1]][[nm2]][[geneNm]] <- list(
            entrez = entrez,
            genesymbol = geneNm,
            weight = 1
          )
          lookup[[nm2]][[nm1]][[geneNm]] <- list(
            entrez = entrez,
            genesymbol = geneNm,
            weight = 1
          )
          # only from large module (src) to small (no need for reverse)
          if(de.num[nm1] >= de.num[nm2]){
            data[[item.count]]<-list(
              fromModule = nm1,
              toModule= nm2,
              arc=geneNm,
              weight=wt,
              fromColor=as.character(colors[nm1]),
              toColor=as.character(colors[nm2])
            );
          }else{
            data[[item.count]]<-list(
              fromModule = nm2,
              toModule=nm1,
              arc=geneNm,
              weight=wt,
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
    weights = weights,
    data = data,
    colors = colors,
    labels = labels
  )
  
  return (
    list(
      chordData = chordData,
      lookup = lookup
    )
  );
}
