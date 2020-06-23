
#' Generate vector of status after one spread step
#'
#' This function generates a vector of status after one spread step.  It assumes that the input adjacency matrix is composed of 0s and 1s
#'
#' Updated 2020-06-03

#' @param amat the adjacency matrix, composed of 0s and 1s
#' @param diag1 if T, diagonal is given value 1 (default is T)
#' @param mattype 'fromto' indicates that rows are sources and columns are sinks, 'tofrom' indicates that columns are sources and rows are sinks (default is 'fromto')
#' @param max1 if T, values in vect2 (the output vector) greater than 1 are replaced by 1 (default is T)
#' @param vect1 status of nodes before spread, in terms of presence or absence of bioentity or information about management (a matrix with nrow=1 and ncol=number of nodes)
#' @keywords dispersal
#' @export 

#' @examples
#' x1 <- spreadstep(amat=matrix(c(1,1,0,0,1,1,0,0,1),ncol=3,byrow=T), vect1=matrix(c(1,0,0), ncol=3), mattype='fromto')
#' x2 <- spreadstep(amat=matrix(c(1,1,0,0,1,1,0,0,1),ncol=3,byrow=T), vect1=matrix(c(1,0,0), ncol=3), mattype='tofrom')
#' x3 <- spreadstep(amat=matrix(c(1,1,0,0,1,1,0,0,1),ncol=3,byrow=T), vect1=matrix(c(1,1,0), ncol=3), mattype='fromto', max1=F)
#' x9 <- spreadstep(amat=matrix(c(0,1,0,0,0,0, 0,0,1,0,0,0, 0,0,0,1,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1, 1,0,0,0,0,0),byrow=T,ncol=6), vect1=c(1,1,1,0,0,0), mattype='fromto')


spreadstep <- function(amat, vect1, diag1=T, mattype='fromto', max1=T) {

  if(diag1) {diag(amat) <- 1}

  if(mattype == 'tofrom'){amat <- t(amat)}
  
  vect2 <- vect1 %*% amat

  if(max1) {vect2 <- as.numeric(vect2 >= 1)}
    
  vect2 # the resulting status of nodes 
}