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
#' @rdname PrepareChordData
#' @export 
PrepareChordData <-function(){
  if(anal.type == "metadata"){
    res <- PrepareMetaChordData();
  }else{
    res <- PrepareListChordData();
  }
  
  if(is.null(res)){
    return(0);
  }
  require(RJSONIO);
  chordData <- res$chordData;
  fileNm = paste0("networkanalyst_chorddata_",chord_count,".json")
  sink(fileNm);
  cat(toJSON(chordData));
  sink();
  
  lookup <- res$lookup;
  fileNm2 = paste0("networkanalyst_chord_lookup_",chord_count,".json")
  sink(fileNm2);
  cat(toJSON(lookup));
  chord_count <<- chord_count + 1
  sink();
  return(1);
}
