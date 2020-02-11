#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param vec PARAM_DESCRIPTION, Default: y.vec
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[BiocGenerics]{which}},\code{\link[BiocGenerics]{colSums}}
#'  \code{\link[stringi]{stri_isempty}},\code{\link[stringi]{stri_length}},\code{\link[stringi]{stri_numbytes}},\code{\link[stringi]{stri_width}}
#'  \code{\link[data.table]{fifelse}}
#' @rdname CalculateS2N
#' @export 
#' @importFrom BiocGenerics which rowMeans
#' @importFrom stringi length
#' @importFrom data.table ifelse
CalculateS2N <- function(data, vec = y.vec, ...) {
  
  A <- data + 0.00000001
  
  ind1 <- BiocGenerics::which(vec==1) # cases
  n1 <- stringi::length(ind1)    
  M1 <- BiocGenerics::rowMeans(A[,ind1])
  A2 <- A*A    
  S1 <- BiocGenerics::rowMeans(A2[,ind1])   
  S1 <- S1 - M1*M1    
  S1 <- sqrt(abs((n1/(n1-1)) * S1))   
  
  ind2 <- BiocGenerics::which(vec==0) # controls
  n2 <- stringi::length(ind2)
  M2 <- BiocGenerics::rowMeans(A[,ind2])
  S2 <- BiocGenerics::rowMeans(A2[,ind2])   
  S2 <- S2 - M2*M2    
  S2 <- sqrt(abs((n2/(n2-1)) * S2))   
  
  # small sigma "fix" as used in GeneCluster
  S2 <- data.table::ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
  S2 <- data.table::ifelse(S2 == 0, 0.2, S2) 
  S1 <- data.table::ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
  S1 <- data.table::ifelse(S1 == 0, 0.2, S1) 
  M1 <- M1 - M2
  S1 <- S1 + S2
  s2n <- M1/S1
  
  return(s2n)
}
