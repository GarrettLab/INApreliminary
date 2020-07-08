
#' Decisions at each node about whether to adopt management
#'
#' This function generates decisions at each node about management
#'
#' Updated 2020-07-07

#' @param adpm mean probability of adopting management if informed
#' @param adpsd sd for drawing probability of adoption, in truncated normal distribution
#' @param comvec vector of 1=info is present, 0=info is not present
#' @param padopt vector of probabilities of adoption if informed for nodes in the network (\code{readpadopt} determines whether \code{padopt} is read in to \code{makedec} or generated within \code{makedec} based on \code{adpm} and \code{adpsd})
#' @param plotmp if T, then map of decision is plotted
#' @param readpadopt if T, then a VECTOR of padopt values is read in; if F, \code{adpm} and \code{adpsd} are used to generate the vector
#' @param xymat matrix of x,y coordinates for the nodes, for mapping if \code{plotmp} = T (and for determining the number of nodes in current version)
#' @keywords decisions
#' @export
#' @import truncnorm

#' @examples
#' x5 <- makedec(comvec=c(1,1,1,0,0,0), xymat=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T,ncol=2), readpadopt=F, adpm=0.1, adpsd=0.1, plotmp=T)
#' x6 <- makedec(comvec=c(1,1,1,0,0,0), xymat=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T,ncol=2), readpadopt=F, adpm=0.5, adpsd=0.1, plotmp=T)
#' x7 <- makedec(comvec=c(1,1,1,0,0,0), xymat=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T,ncol=2), readpadopt=T, padopt=c(0.1,0.1,0.9,0.1,0.2,0.3), plotmp=T)


makedec <- function(comvec, xymat, adpm, adpsd, readpadopt, padopt, plotmp=F){

  nodlen <- dim(xymat)[1] # number of nodes
 
  if(!readpadopt){
    padopt <- truncnorm::rtruncnorm(n=nodlen, a=0, b=1, mean = adpm, sd=adpsd)
  }
  
  stochadopt <- runif(n=nodlen, min=0, max=1)
  hypdecvec <- stochadopt < padopt # would adopt if informed
  decvec <- hypdecvec & comvec # does adopt

  if(plotmp){
    plot(xymat, main='Decision to adopt, shaded = yes')
    if(!is.null(dim(xymat[decvec==1,]))){
      points(xymat[decvec==1,], pch=16)
    } 
    else if(is.null(dim(xymat[decvec==1,]))){
      points(xymat[decvec==1,][1], xymat[decvec==1,][2], pch=16)
    }
  }

  decvec
}