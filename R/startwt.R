

#' Evaluate nodes for value in invasion detection (based on a specified number of realizations, and including specified probabilities associated with each potential starting node)
#'
#' Uses weights that indicate how likely each node was to be the starting location (introduction point) for an invasion.  Uses output from multisimtest.  Finds the rate at which sampling success (other nodes free of invasion when invasion reaches sampling node) is expected to occur for each potential sampling node.
#' @param msf.out output object from function multisimf
#' @param adjmat adjacency matrix for evaluation
#' @param wtvec vector of weights indicating the probability that each node would be the starting node for invasion
#' @param nodenam vector of node names
#' @keywords prioritization sampling weights
#' @export
#' @examples
#' wtvec.ex <- c(1,10,100,1)
#' msf.outex <- multisimtest(adjmat=sAmat, stoch=T, nsim=10)
#' startwt(msf.out=msf.outex, adjmat=sAmat, wtvec=wtvec.ex)
#' startwt(msf.out=msf.outex, adjmat=sAmat, wtvec=wtvec.ex, nodenam=c("KS","NE","ND","SD"))
#' startwt(msf.out=msf.outex, adjmat=sAmat, wtvec=c(1,100,10,1), nodenam=c("KS","NE","ND","SD"))

# to do - improve examples


startwt <- function(msf.out, adjmat, wtvec, nodenam=NA) {

  dimL <- dim(adjmat)[1]

  matop <- msf.out$meanarr
  wtarr <- wtvec * matop

  # find invasion free rate for each sampling node
  sampfree <- colSums(wtarr)
  tsampfree <- data.frame(1:dimL, sampfree, nodenam)
  
  list(wtarr=wtarr, tsampfree=tsampfree)

}

