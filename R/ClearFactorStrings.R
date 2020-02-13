#' @title FUNCTION_TITLE
#' @description utils to remove from within
#' @param cls.nm PARAM_DESCRIPTION
#' @param query PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details removes leading and trailing spaces
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname ClearFactorStrings
#' @export 
ClearFactorStrings<-function(cls.nm, query){
  # remove leading and trailing space
  query<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", query, perl=TRUE);
  
  # kill multiple white space
  query <- gsub(" +","_",query);
  # remove non alphabets and non numbers 
  query <- gsub("[^[:alnum:] ]", "_", query);
  
  # test all numbers (i.e. Time points)
  chars <- substr(query, 0, 1);
  num.inx<- chars >= '0' & chars <= '9';
  if(all(num.inx)){
    query = as.numeric(query);
    nquery <- paste(cls.nm, query, sep="_");
    query <- factor(nquery, levels=paste(cls.nm, sort(unique(query)), sep="_"));
  }else{
    query[num.inx] <- paste(cls.nm, query[num.inx], sep="_");
    query <- factor(query);
  }
  return (query);
}
