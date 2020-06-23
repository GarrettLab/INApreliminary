
#' Generate network adjacency matrix for movement
#'
#' This function generates an adjacency matrix for movement, assumed symmetric (in this version).  It is used by functions including \code{INAscene}. The movement adjacency matrix is composed of 1s and 0s only if lktype="pa" option is used
#'
#' Updated 2020-06-02

#' @param distf the function of distance used to estimate movement probability - 'random' (not related to distance) or 'powerlaw' (inverse power law) or 'exp' (negative exponential, to be added)
#' @param iplot if T, generates igraph plot of adjacency matrix
#' @param lktype link type, pa is presence/absence (unweighted, occurence/non-occurence), pr is a probability of occurence, wtd1 is general weight
#' @param pla inverse power law parameter a in ad^(-b)
#' @param plb inverse power law parameter b in ad^(-b)
#' @param randp random matrix with entries binomial with probability p
#' @param tlink threshold for whether link exists 
#' @param xymat the matrix of xy coordinates for node locations, used when the probability of a link is a function of distance (note that the distance between each pair of locations is assumed to be greater than 1)
#' @keywords dispersal
#' @export 

#' @examples
#' x1 <- genmovnet(j <- genlocs(extx=c(0,50), exty=c(0,50), nn=50, rand=TRUE), distf='random', randp=0.01, lktype='pa', tlink=0.05, iplot=T)
#' x2 <- genmovnet(j <- genlocs(extx=c(0,50), exty=c(0,50), nn=100, rand=TRUE), distf='random', randp=0.02, lktype='pa', tlink=0.05, iplot=T)
#' x7 <- genmovnet(xymat=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3),ncol=2,byrow=T), distf='powerlaw', pla=2, plb=1, lktype='pa', tlink=0.9, iplot=T)
#' x8 <- genmovnet(j <- genlocs(nn=30, extx = c(0, 10), exty = c(0, 10)), distf='powerlaw', pla=2, plb=1, lktype='pa', tlink=0.9, iplot=T)
#' x9 <- genmovnet(j <- genlocs(nn=300, extx = c(0, 10), exty = c(0, 100)), distf='powerlaw', pla=2, plb=1, lktype='pa', tlink=0.95, iplot=T)


genmovnet <- function(xymat, distf, iplot=F, lktype, randp, pla, plb, tlink){

  dimam <- dim(xymat)[1]

  if (distf == 'powerlaw') { # ad^(-b)

    tdist <- as.matrix(dist(xymat, method = "euclidean", diag=T, upper=T))
    linkmat <- pla*tdist^(-plb)  
  }

  else if(distf == 'random'){
  
    linkmat <- matrix(rbinom(n=dimam*dimam, size=1, prob=randp), nrow=dimam)

    linkmat[lower.tri(linkmat)] <- t(linkmat)[lower.tri(linkmat)]
    
  }

  
# If generating presence/absence of links for cases other than 'random', keep links where linkmat entries are above the threshold tlink

  if(lktype == 'pa' & distf == 'powerlaw'){linkmat <- linkmat > tlink}

  diag(linkmat) <- 0  


  if(iplot){ 
    
    library(igraph)
    linkmati <- graph.adjacency(linkmat)
    plot(linkmati, edge.arrow.size=0, vertex.color='skyblue')
    
  }

  linkmat
}