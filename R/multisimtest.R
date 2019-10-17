  
#' Evaluate nodes for value in invasion detection (based on a specified number of realizations)
#'
#' Evaluate nodes for utility in invasion detection (based on a specified number of realizations) using the function multistart (which in turn uses the function onestart).  outarr saves the output from each realization from multistart, in a 3D array.  In meanarr, rows = starting nodes (introduction nodes), columns = sampling nodes, and entries = mean nodes free from invasion when invasion starting at the introduction node reaches the sampling node.  In vararr, for the same rows and columns, entries are the variance in the number of nodes free from invasion.  Note that there is no reason to use multisimtest if there is no stochastic component.
#' @param adjmat adjacency matrix for evaluation
#' @param stoch logical var indicating whether adjacency matrix entries are fixed or probabilities
#' @param nsim number of simulations to be analyzed
#' @keywords prioritization sampling
#' @export
#' @examples
#' multisimtest(adjmat=sAmat, stoch=T, nsim=10)[

# to do - consider skewedness?  kurtosis?
# to do - note that can use stoch=F, for output to startwt (or change startwt to also use output from multistart)


multisimtest <- function(adjmat, stoch, nrealz=1) {

  dimL <- dim(adjmat)[1]

  outarr <- array(-99, c(dimL,dimL,nrealz)) # all results
  meanarr <- matrix(-99, ncol=dimL, nrow=dimL) # mean of results
  vararr <- meanarr # variance of results

  for (i3 in 1:nrealz){
    outarr[,,i3]  <- multistart(adjmat=adjmat, stoch=stoch)
  }

  for (i3 in 1:dimL) {
    for (i4 in 1:dimL) {
      meanarr[i3,i4] <- mean(outarr[i3,i4,]) #mean across realizations
      vararr[i3,i4] <- var(outarr[i3,i4,]) #variance across realizations
  } }
# to do - consider alternative to replace the second loop

  list(outarr=outarr, meanarr=meanarr, vararr=vararr)

}


