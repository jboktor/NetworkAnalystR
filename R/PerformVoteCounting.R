# diff used for direction, not selection
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param BHth PARAM_DESCRIPTION, Default: 0.05
#' @param minVote PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname PerformVoteCounting
#' @export 
PerformVoteCounting <- function(BHth = 0.05, minVote){
  if(!performedDE){
    PerformMetaDeAnal();
  }
  inmex.meta <- readRDS("inmex_meta.rds");
  inmex.method <<- "votecount";
  DE.vec <<- NULL; # store entrez id from meta-analysis for GO
  meta.mat <<- meta.stat <<- NULL;
  sel.nms <- names(mdata.all)[mdata.all==1];
  # first create a matrix to stall the result
  # row for each feature and col for each dataset uploaded
  vc.mat <- matrix(0, nrow=nrow(inmex.meta$data), ncol=length(sel.nms)+1);
  shared.ids <- rownames(inmex.meta$data);
  for(i in 1:length(inmex.ind)){
    res.mat <- inmex.ind[[i]];
    res.mat <- res.mat[shared.ids, ];
    
    #note in meta-analysis should consider directions
    # use logFC for this purpose 
    # consider upregulated
    hit.up.inx <- res.mat[,1]> 0 & res.mat[,2] <= BHth;
    up.vote <- as.numeric(hit.up.inx);
    
    # consider downregulated
    hit.dn.inx <- res.mat[,1] < 0 & res.mat[,2] <= BHth;
    dn.vote <- -as.numeric(hit.dn.inx);
    
    vc.mat[,i] <- up.vote + dn.vote;
  }
  
  # total score (votes for each direction)
  vc.mat[,length(sel.nms)+1] <- apply(vc.mat, 1, sum);
  colnames(vc.mat) <- c(paste("Vote", substring(sel.nms,0, nchar(sel.nms)-4)), "VoteCounts");
  rownames(vc.mat) <- rownames(inmex.meta$data);
  
  # compute at least one vote (no direction)
  vote.any <- apply(abs(vc.mat), 1, sum)
  vote.any.inx <- vote.any > 0;
  
  # return results with at least one vote
  vc.mat <- vc.mat[vote.any.inx, ];
  
  #sort
  ord.inx <- order(abs(vc.mat[, "VoteCounts"]), decreasing = T);
  vc.mat <- vc.mat[ord.inx, "VoteCounts", drop=F];
  
  sig.inx <- abs(vc.mat[,"VoteCounts"]) >= minVote;
  meta.mat <<- vc.mat;
  meta.mat.all <<- vc.mat;
  SetupMetaStats(BHth);
  return(length(sig.inx));
}
