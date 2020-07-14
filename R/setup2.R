
#' Set up starting conditions for INA scenario analysis
#'
#' This function sets up all starting conditions, including the estimated effect size, the x,y coordinates of the nodes, the initial information locations, and the initial species locations. It runs functions \code{estinfo}, \code{genlocs}, and \code{initvals}. The next step after this is to address the adjacency matrices for communication and for dispersal.
#'
#' Updated 2020-07-12

#' @param manredestab2 if T, the management reduces the probability of establishment by probability with mean efmn2 and sd efsd2; if F, the management reduces the probability of NOT establishing by the same
#' @param efmn2 (estinfo) the underlying mean change in establishment probability (as a proportion)
#' @param efsd2 (estinfo) the standard deviation of the effect
#' @param efthr2 (estinfo) the threshold effect size for communicating about management (if efthr2 = 0 there is no threshold so communication can always occur)
#' @param sampeffort2 (estinfo) sampling effort, where greater samping effort reduces the error in estimating the management effect

#' @param readxymatj read in the xy coordinates for locations (readxymatj = T) or generate the xy coordinates by calling genlocs (readxymatj = F)
#' @param xymat2 coordinates for node geographic locations, ncol= 2 (x,y) and nrow=number of nodes (to be read in if readorgenxy = 'read')

#' @param extx2 (genlocs) range of x coordinates
#' @param exty2 (genlocs) range of y coordinates
#' @param nodnum2 (genlocs) the number of nodes (nn in genlocs)
#' @param rand2 (genlocs) if TRUE then locations are randomly generated

#' @param readinitinfo.s if T, the initial values for the vector of starting locations for the presence of information are read in rather than generated
#' @param initinfo.s the vector of initial values read in if readinitinfo.s == T, a vector with lenth equal to the number of nodes and entries 1s or 0s with 1s indicating the initial presence of information

#' @param loctypei info (initvals) the type of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence
#' @param numiniti info (initvals) the number of initial locations for presence
#' @param numorpropi info (initvals) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param propiniti info (initvals) the proportion of initial locations for presence (ip in initvals)

#' @param readinitbio.s if T, the initial values for the vector of starting locations for the presence of the bioentity are read in rather than generated
#' @param initbio.s the vector of initial values read in if readinitbio.s == T, a vector with lenth equal to the number of nodes and entries 1s or 0s with 1s indicating the initial presence of the bioentity

#' @param loctypeb bioentity (initvals) the type of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence
#' @param numinitb bioentity (initvals) the number of initial locations for presence
#' @param numorpropb bioentity (initvals) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param propinitb bioentity (initvals) the proportion of initial locations for presence (ip in initvals)

#' @param plotmp if T plots of initial presence of information and species are generated

#' @keywords simulation experiments
#' @export 

#' @examples
#' x2 <- setup2(manredestab2=T, efmn2=0.5, efsd2=0.5, efthr2=0.5, sampeffort2=1, usethresh2=T, readxymatj=F, extx2=c(0,50), exty2=c(0,50), nodnum2=100, rand2=TRUE, readinitinfo.s=F, loctypei='random', numiniti=5, numorpropi='num', readinitbio.s=F, loctypeb='upedge', numinitb=5, numorpropb='num', plotmp=T)
#' x3 <- setup2(manredestab2=T, efmn2=0.5, efsd2=0.5, efthr2=0.5, sampeffort2=1, usethresh2=T, readxymatj=F, xtx2=c(0,50), exty2=c(0,50), nodnum2=100, rand2=TRUE, readinitinfo.s=F, loctypei='random', numiniti=15, numorpropi='num', readinitbio.s=F, loctypeb='upedge', numinitb=15, numorpropb='num', plotmp=T)
#' x4.readxy <- setup2(manredestab2=T, efmn2=0.5, efsd2=0.5, efthr2=0.5, sampeffort2=1, usethresh2=T, readxymatj=T, xymat2=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3),byrow=T,ncol=2), readinitinfo.s=F, loctypei='random', numiniti=2, numorpropi='num', readinitbio.s=F, loctypeb='upedge', numinitb=3, numorpropb='num', plotmp=T)


setup2 <- function(manredestab2, efmn2, efsd2, efthr2, sampeffort2, usethresh2, readxymatj, xymat2, extx2, exty2, nodnum2, rand2, readinitinfo.s, initinfo.s, loctypei, numiniti, numorpropi, propiniti, readinitbio.s, initbio.s, loctypeb, numinitb, numorpropb, propinitb, plotmp=F){

if (usethresh2) {
  infout <- estinfo(efmn=efmn2, efsd=efsd2, efthr=efthr2,   sampeffort=sampeffort2)
} else {
  infout <- list()
  infout$com.yes <- T
  infout$obschange <- NA
}


if (readxymatj == F){
  xymat2 <- genlocs(extx=extx2, exty=exty2, nn=nodnum2, rand=rand2)
}

# generate initial locations for info if not supplied
if (readinitinfo.s == FALSE) {
  infovec <- initvals(xymat=xymat2, loctype=loctypei,    numinit=numiniti, numorprop=numorpropi, ip=propiniti) 
} else if (readinitinfo.s == TRUE) {
  infovec <- initinfo.s
}

# generate initial locations for bioentity establishment
#   if not supplied
if (readinitbio.s == FALSE) {
  estabvec <- initvals(xymat=xymat2, loctype=loctypeb, numinit=numinitb, numorprop=numorpropb, ip=propinitb) 
} else if (readinitbio.s == TRUE) {
  estabvec <- initbio.s
}

if(plotmp){
    plot(xymat2, main='Information, shaded = yes')
    if(!is.null(dim(xymat2[infovec==1,]))){
      points(xymat2[infovec==1,], pch=16)
    } 
    else if(is.null(dim(xymat2[infovec==1,]))){
      points(xymat2[infovec==1,][1], xymat2[infovec==1,][2], pch=16)
    }
  }

if(plotmp){
    plot(xymat2, main='Bioentity establishment, shaded = yes')
    if(!is.null(dim(xymat2[estabvec==1,]))){
      points(xymat2[estabvec==1,], pch=16)
    } 
    else if(is.null(dim(xymat2[estabvec==1,]))){
      points(xymat2[estabvec==1,][1], xymat2[estabvec==1,][2], pch=16)
    }
  }

  list(com.yes=infout$com.yes, obschange=infout$obschange, xymat2=xymat2, infovec=infovec, estabvec=estabvec)
}