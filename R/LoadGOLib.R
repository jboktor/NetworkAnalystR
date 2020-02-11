#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param onto PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname LoadGOLib
#' @export 
LoadGOLib<-function(onto){
  go.path <- paste(lib.path, data.org, "/go_", tolower(onto), ".rds", sep="");
  if(tolower(onto) == "bp"){
    go_bp <- readRDS(go.path);
    if(is.null(names(go_bp))){ # new go lib does not give names
      names(go_bp) <- c("link", "term", "sets");
    }
    current.link <- go_bp$link;
    current.geneset <- go_bp$sets;
    set.ids<- names(current.geneset); 
    names(set.ids) <- names(current.geneset) <- go_bp$term;
  }else if(tolower(onto) == "mf"){
    go_mf <- readRDS(go.path);
    if(is.null(names(go_mf))){
      names(go_mf) <- c("link", "term", "sets");
    }
    current.link <- go_mf$link;
    current.geneset <- go_mf$sets;
    set.ids<- names(current.geneset); 
    names(set.ids) <- names(current.geneset) <- go_mf$term;
  }else{
    go_cc <- readRDS(go.path);
    if(is.null(names(go_cc))){
      names(go_cc) <- c("link", "term", "sets");
    }
    current.link <- go_cc$link;
    current.geneset <- go_cc$sets;
    set.ids<- names(current.geneset); 
    names(set.ids) <- names(current.geneset) <- go_cc$term;
  }
  names(current.geneset) = firstup(names(current.geneset))
  names(current.geneset) = gsub("-", "_", names(current.geneset))
  names(set.ids) = firstup(names(set.ids));
  names(set.ids) = gsub("-", "_", names(set.ids))
  current.setlink <<- current.link;
  current.setids <<- set.ids;
  saveRDS(current.geneset, "current_geneset.rds");
  return(current.geneset);
}
