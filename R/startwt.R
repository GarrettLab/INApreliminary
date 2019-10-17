

#' Evaluate nodes for value in invasion detection (based on a specified number of simulations, and including specified probabilities associated with each potential starting node)
#'
#' Uses weights that indicate how likely each node was to be the starting location for an invasion.  Uses output from multisimtest.  Finds the rate at which sampling failure (invasion of other nodes prior to detection at sampling node) is expected to occur for each potential sampling node.
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

# to do - GT
# to do - improve examples
#

startwt <- function(msf.out, adjmat, wtvec, nodenam=NA) {

  dimL <- dim(adjmat)[1]

  matop <- msf.out$meanarr
  wtarr <- wtvec * matop

  # find failure rate for each sampling node
  sampfail <- colSums(wtarr)
  tsampfail <- data.frame(1:dimL, sampfail, nodenam)
  
  list(wtarr=wtarr, tsampfail=tsampfail)

}

