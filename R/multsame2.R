
#' Runs and summarizes multiple INA simulations with the same parameter values
#'
#' This function runs and summarizes multiple INA simulations with the same parameter values.  (It uses functions estinfo, genlocs, initvals, setup2, genmovnet, spreadstep, makedec, estab, and ntsteps2.)
#'
#' Updated 2020-06-11

#' @param nreals number of realizations to be evaluated
#' @param usethresh if T, information is never present anywhere unless the management effect estimate exceeds the threshold
#' @param nsteps2 number of time steps to be evaluated

#' @param readcomamn2 if T, the communication adjacency matrix is read in rather than being generated from the outset
#' @param comamn2 communication adjacency matrix, read in if readcomamn2=T
#' @param readdispamn2 if T, the dispersal adjacency matrix is read in rather than being generated from the outset
#' @param dispamn2 dispersal adjacency matrix, read in if readdispamn2=T

#' @param efmn22 (estinfo) the underlying mean change in establishment probability (as a proportion)
#' @param efsd22 (estinfo) the standard deviation of the effect
#' @param efthr22 (estinfo) the threshold effect size for communicating about management (if efthr22 = 0 there is no threshold so communication can always occur)
#' @param sampeffort22 (estinfo) sampling effort, where greater samping effort reduces the error in estimating the management effect

#' @param readorgenxy22 read in the xy coordinates for locations (readorgenxy22 = 'read') or generate the xy coordinates by calling genlocs (readorgenxy22 = 'gen')

#' @param readxymat2 if T, read in xymatinn2 - otherwise, generate it in each realization (related readxymat2=T and readorgenxy='read')
#' @param xymatinn2 matrix of x,y coordinates of nodes, read in if readxymat2=T

#' @param extx22 (genlocs) range of x coordinates
#' @param exty22 (genlocs) range of y coordinates
#' @param nodenum22 (genlocs) the number of nodes
#' @param rand22 (genlocs) if TRUE then locations are randomly generated

#' @param readinitinfo.s2 info (setup2) if T, the initial values for the vector of starting locations for the presence of information are read in rather than generated
#' @param initinfo.s2 info (setup2) the vector of initial values read in if readinitinfo.s == T, a vector with lenth equal to the number of nodes and entries 1s or 0s with 1s indicating the initial presence of information

#' @param loctypei2 info (initvals) the type of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence
#' @param numiniti2 info (initvals) the number of initial locations for presence
#' @param numorpropi2 info (initvals) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param propiniti2 info (initvals) the proportion of initial locations for presence

#' @param readinitbio.s2 bio (setup2) if T, the initial values for the vector of starting locations for the presence of the bioentity are read in rather than generated
#' @param initbio.s2 bio (setup2) the vector of initial values read in if readinitbio.s == T, a vector with lenth equal to the number of nodes and entries 1s or 0s with 1s indicating the initial presence of the bioentity

#' @param loctypeb2 estab (initvals) the type of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence
#' @param numinitb2 estab (initvals) the number of initial locations for presence
#' @param numorpropb2 estab (initvals) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param propinitb2 estab (initvals) the proportion of initial locations for presence

#### ending section of startup components

#' @param distfcn2 com (genmovnet) the function of distance used to estimate movement probability - 'random' (not related to distance) or 'powerlaw' 
#' @param lktypecn2 com (genmovnet) link type, pa is presence/absence (unweighted, occurence/non-occurence), pr is a probability of occurence, wtd1 is general weight
#' @param placn2 com (genmovnet) power law parameter a in ad^(-b)
#' @param plbcn2 com (genmovnet) power law parameter b in ad^(-b)
#' @param randpcn2 com (genmovnet) random case, probability of link
#' @param tlinkcn2 com (genmovnet) threshold for whether link exists (probability in some cases); used if lktypecn2 = pa

#' @param distfdn2 disp (genmovnet) the function of distance used to estimate movement probability - 'random' (not related to distance) or 'powerlaw' 
#' @param lktypedn2 disp (genmovnet) link type, pa is presence/absence (unweighted, occurence/non-occurence), pr is a probability of occurence, wtd1 is general weight
#' @param pladn2 disp (genmovnet) power law parameter a in ad^(-b)
#' @param plbdn2 disp (genmovnet) power law parameter b in ad^(-b)
#' @param randpdn2 disp (genmovnet) random case, probably of link
#' @param tlinkdn2 disp (genmovnet) threshold for whether link exists (probability in some cases); used if lktypedn2 = pa (to be confirmed)

#' @param diag1cn2 com (spreadstep) if T, diagonal of adj mat is given value 1
#' @param mattypecn2 com (spreadstep) 'fromto' indicates that rows are sources and columns are sinks, 'tofrom' indicates that columns are sources and rows are sinks
#' @param max1cn2 com (spreadstep) if T, values in vect2 greater than 1 are replaced by 1
#' @param vect1cn2 com (spreadstep) status of nodes before spread

#' @param diag1dn2 disp (spreadstep) if T, diagonal of adj mat is given value 1
#' @param mattypedn2 disp (spreadstep) 'fromto' indicates that rows are sources and columns are sinks, 'tofrom' indicates that columns are sources and rows are sinks
#' @param max1dn2 disp (spreadstep) if T, values in vect2 greater than 1 are replaced by 1
#' @param vect1dn2 disp (spreadstep) status of nodes before spread

#' @param adpmn2 (makedec) mean probability of adopting management if informed
#' @param adpsdn2 (makedec) sd in truncated normal distribution of probability of adoption
#' @param padopt2 vector of probabilities of adoption if informed for nodes in the network (\code{readpadopt} determines whether \code{padopt} is read in to \code{makedec} or generated within \code{makedec} based on \code{adpm} and \code{adpsd})
#' @param readpadopt2 if T, then a VECTOR of padopt values is read in; if F, \code{adpm} and \code{adpsd} are used to generate the vector

#' @param estpmn2 (estab) mean probability of establishment (new or CONTINUED) in absence of management
#' @param estpsdn2 (estab) sd in truncated normal distribution for probability of establishment in absence of management
#' @param efmnn2 (estab) mean management effect (proportional reduction in estabp) (same as efmn22 to function estinfo)
#' @param efsdn2 (estab) sd of management effect (same as efsd22 to function estinfo)
#' @param pestabn2 (estab) vector of probabilities of establishment (read in or generated when the function is run, depending on \code{readpestab} equal to T or F)
#' @param readpestabn2 (estab) if T, then a vector \code{pestabn2} is read in; otherwise the vector is generated using \code{estpmn2} and \code{estpsdn2}

#' @param plotmpn2 if true plots of resulting presence of information and species are generated

#' @keywords simulation experiments
#' @export
#' @examples


#' x5 <- multsame2(nreals=10, nsteps2=3, usethresh=F, readxymat2=T, readorgenxy22='read', xymatinn2=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3),byrow=T,ncol=2), readcomamn2=F, readdispamn2=F, efmn22=0.5, efsd22=0.1, efthr22=0.5, sampeffort22=1, extx22=NA, exty22=NA, nodenum22=NA, rand22=NA, readinitinfo.s2=F, loctypei2='random', numiniti2=5, numorpropi2='num', propiniti2=NA, readinitbio.s2=F, loctypeb2='upedge', numinitb2=5, numorpropb2='num', propinitb2=NA, distfcn2='powerlaw', lktypecn2='pa', placn2=1, plbcn2=1, randpcn2=NA, tlinkcn2=0.1, distfdn2='powerlaw', lktypedn2='pa', pladn2=1, plbdn2=1, randpdn2=NA, tlinkdn2=0.1, diag1cn2=T, mattypecn2='fromto', max1cn2=T, diag1dn2=T, mattypedn2='fromto', max1dn2=T, adpmn2=0.5, adpsdn2=0.1, padopt2=NA, readpadopt2=F, estpmn2=0.5, estpsdn2=0.1, pestabn2=NA, readpestabn2=F, plotmpn2=F)

#' x3 <- multsame2(nreals=2,nsteps=3, usethresh=F, readxymat2=T, readorgenxy22='read', xymatinn2=matrix(runif(n=100)*100,byrow=T,ncol=2), readcomamn2=F, readdispamn2=F, efmn22=0.5, efsd22=0.1, efthr22=0.5, sampeffort22=1, extx22=NA, exty22=NA, nodenum22=NA, rand22=NA, readinitinfo.s2=F, loctypei2='random', numiniti2=5, numorpropi2='num', propiniti2=0.05, readinitbio.s2=F, loctypeb2='upedge', numinitb2=5, numorpropb2='num', propinitb2=0.05, distfcn2='powerlaw', lktypecn2='pa', placn2=1, plbcn2=1, randpcn2=NA, tlinkcn2=0.1, distfdn2='powerlaw', lktypedn2='pa', pladn2=1, plbdn2=1, randpdn2=NA, tlinkdn2=0.1, diag1cn2=T, mattypecn2='fromto', max1cn2=T, diag1dn2=T, mattypedn2='fromto', max1dn2=T, adpmn2=0.5, adpsdn2=c(0,1), padopt2=NA, readpadopt2=F, estpmn2=0.5, estpsdn2=0.1, pestabn2=NA, readpestabn2=F, plotmpn2=F)



multsame2 <- function(nreals, usethresh, nsteps2, xymatinn2, readxymat2, readorgenxy22, comamn2, dispamn2, readcomamn2, readdispamn2, efmn22, efsd22, efthr22, sampeffort22, extx22, exty22, nodenum22, rand22, readinitinfo.s2, initinfo.s2, loctypei2, numiniti2, numorpropi2, propiniti2, readinitbio.s2, initbio.s2, loctypeb2, numinitb2, numorpropb2, propinitb2, distfcn2, lktypecn2, placn2, plbcn2, randpcn2, tlinkcn2, distfdn2, lktypedn2, pladn2, plbdn2, randpdn2, tlinkdn2, diag1cn2, mattypecn2, max1cn2, diag1dn2, mattypedn2, max1dn2, adpmn2, adpsdn2, padopt2, readpadopt2, estpmn2, estpsdn2, pestabn2, readpestabn2, plotmpn2=F){



  if(readxymat2){nodenum22 <- dim(xymatinn2)[1]}

  # prepare an output object with output from ntsteps2 
  #   for each realization
  multout <- as.list(1:nreals)
  
  # prepare an output object with output from setup2 
  #   for each realization
  setupout <- multout

  # keep the proportion of nodes in each category:
  #   with com, with adoption, with disp, with establishment
  setocom <- 0*(1:nreals)
  setodec <- setocom
  setodisp <- setocom
  setoestab <- setocom

  for(j in 1:nreals) {


    tempo <- setup2(efmn2=efmn22, efsd2=efsd22, efthr2=efthr22, sampeffort2=sampeffort22, usethresh2=usethresh, readorgenxy=readorgenxy22, xymat2=xymatinn2, extx2=extx22, exty2=exty22, nodnum2=nodenum22, rand2=rand22, readinitinfo.s=readinitinfo.s2, initinfo.s=initinfo.s2, loctypei=loctypei2, numiniti=numiniti2, numorpropi=numorpropi2, propiniti=propiniti2, readinitbio.s=readinitbio.s2, initbio.s=initbio.s2, loctypeb=loctypeb2, numinitb=numinitb2, numorpropb=numorpropb2, propinitb=propinitb2, plotmp=plotmpn2)

# note that the following draws on tempo

    temp2 <- ntsteps2(nsteps=nsteps2, infon=tempo$com.yes, xymatinn=tempo$xymat2, vect1cn=tempo$infovec, vect1dn=tempo$estabvec, comamn=comamn2, dispamn=dispamn2, readcomamn=readcomamn2, readdispamn=readdispamn2, distfcn=distfcn2, lktypecn=lktypecn2, placn=placn2, plbcn=plbcn2, randpcn=randpcn2, tlinkcn=tlinkcn2, distfdn=distfdn2, lktypedn=lktypedn2, pladn=pladn2, plbdn=plbdn2, randpdn=randpdn2, tlinkdn=tlinkdn2, diag1cn=diag1cn2, mattypecn=mattypecn2, max1cn=max1cn2, diag1dn=diag1dn2, mattypedn=mattypedn2, max1dn=max1dn2, adpmn=adpmn2, adpsdn=adpsdn2, estpmn=estpmn2, estpsdn=estpsdn2, efmnn=efmn22, efsdn=efsd22, pestabn=pestabn2, readpestabn=readpestabn2, plotmpn=plotmpn2)

    # save the output from setup2 and ntsteps2
    multout[[j]] <- temp2
    setupout[[j]] <- tempo

    # the proportion of nodes in each category
    setocom[j] <- sum(temp2$vect1cL[[nsteps2+1]])/nodenum22
    setodec[j] <- sum(temp2$decvecL[[nsteps2]])/nodenum22
    setodisp[j] <- sum(temp2$vect1dL[[nsteps2+1]])/nodenum22
    setoestab[j] <- sum(temp2$estabvecL[[nsteps2]])/nodenum22

  } ### end of nreals realizations

  meancom <- mean(setocom)
  meandec <- mean(setodec)
  meandisp <- mean(setodisp)
  meanestab <- mean(setoestab)

# add the 5th percentile and 95th percentile

  com5 <- quantile(setocom, probs=0.05)
  dec5 <- quantile(setodec, probs=0.05)
  disp5 <- quantile(setodisp, probs=0.05)
  estab5 <- quantile(setoestab, probs=0.05)

  com95 <- quantile(setocom, probs=0.95)
  dec95 <- quantile(setodec, probs=0.95)
  disp95 <- quantile(setodisp, probs=0.95)
  estab95 <- quantile(setoestab, probs=0.95)


  list(multout=multout, setocom=setocom, setodec=setodec, setodisp=setodisp, setoestab=setoestab, meancom=meancom, meandec=meandec, meandisp=meandisp, meanestab=meanestab, com5=com5, dec5=dec5, disp5=disp5, estab5=estab5, com95=com95, dec95=dec95, disp95=disp95, estab95=estab95)
}