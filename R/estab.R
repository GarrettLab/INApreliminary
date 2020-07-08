
#' Bioentity establishment or not at each node
#'
#' This function determines whether a bioentity establishes or not at each node.  It is similar in some ways to the function \code{makedec}, which determines whether management is adopted.
#'
#' Updated 2020-07-07

#' @param decvec vector of 1=management, 0=no management
#' @param dispvec vector of 1=sp has dispersed to node, 0=sp has not dispersed to node
#' @param estpm mean probability of establishment (new or CONTINUED) in absence of management (used if \code{readpestab} =  F)
#' @param estpsd sd in probability of establishment in absence of management (used if \code{readpestab} = F)
#' @param efmn mean effect of management (proportional reduction in estabp, where the combination of efmn = 1 and efsd = 0 makes establishment impossible)
#' @param efsd sd of management effect
#' @param pestab vector of probabilities of establishment (read in or generated when the function is run, depending on \code{readpestab} equal to T or F)
#' @param plotmp if T, then map of establishment is plotted
#' @param readpestab if T, then a vector \code{pestab} is read in; otherwise the vector is generated using \code{estpm} and \code{estpsd}
#' @param xymat matrix of x,y coordinates of nodes (used if \code{plotmp} = T)
#' @keywords establishment
#' @export
#' @import truncnorm

#' @examples
#' x6 <- estab(decvec=c(1,1,1,0,0,0), dispvec=c(1,1,1,1,1,0), xymat=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T, ncol=2), readpestab=F, estpm=0.5, estpsd=0.1, efmn=0.1, efsd=0.1, plotmp=T)
#' x7 <- estab(decvec=c(1,1,1,0,0,0), dispvec=c(1,1,1,1,1,0), xymat=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T, ncol=2), readpestab=F, estpm=0.5, estpsd=0.1, efmn=0.9, efsd=0.1, plotmp=T)
#' x8 <- estab(decvec=c(1,1,1,0,0,0), dispvec=c(1,1,1,1,1,0), xymat=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T, ncol=2), readpestab=T, pestab=c(0.9,0.1,0.1,0.1,0.1,0.9), efmn=0.9, efsd=0.1, plotmp=T)
#' x9 <- estab(decvec=c(1,1,1,0,0,0), dispvec=c(1,1,1,1,1,0), xymat=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3), byrow=T, ncol=2), readpestab=T, pestab=c(0.9,0.1,0.1,0.1,0.1,0.9), efmn=0.1, efsd=0.1, plotmp=T)


estab <- function(decvec, dispvec, estpm, estpsd, efmn, efsd, plotmp=F, xymat, readpestab, pestab){

  nodlen <- length(decvec)

  # probability of establishment in absence of management
  if(!readpestab){
    pestab <- truncnorm::rtruncnorm(n=nodlen, a=0, b=1, mean=estpm, sd=estpsd)
  }

  # how great would the management effect be if used
  obsman <- truncnorm::rtruncnorm(n=nodlen, a=0, b=1, mean=efmn, sd=efsd)

  # probabiliy of estab adjusted for whether management used
  pestab[decvec==1] <- (1-obsman)[decvec==1] * pestab[decvec==1]

  # random draw used for comparison to pestab
  stochestab <- runif(n=nodlen, min=0, max=1)

  # whether would establish if has dispersed
  estabvec <- stochestab < pestab

  # has it also dispersed there?
  estabvec <- estabvec & dispvec

  if(plotmp){
    plot(xymat, main='Establishment, shaded = yes')
    if(!is.null(dim(xymat[estabvec==1,]))){
      points(xymat[estabvec==1,], pch=16)
    } 
    else if(is.null(dim(xymat[estabvec==1,]))){
      points(xymat[estabvec==1,][1], xymat[estabvec==1,][2], pch=16)
    }
  }

  estabvec
}