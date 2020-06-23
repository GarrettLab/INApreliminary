
#' Runs and summarizes multiple simulations with varying parameter values
#'
#' This function runs and summarizes multiple simulations with varying parameter values
#'
#' Updated 2019-08-11

#' @param nsims3 number of simulations to be evaluated
#' @param usethresh3 if T, information is never present anywhere unless the management effect estimate exceeds the threshold
#' @param nsteps3 number of time steps to be evaluated

#' @param readcomamn3 if T, the communication adjacency matrix is read in rather than being generated from the outset
#' @param readdispamn3 if T, the dispersal adjacency matrix is read in rather than being generated from the outset
#' @param xymatn3 matrix of x,y coordinates of nodes

#' @param readxymat3 if T, read in xymatinn2 - otherwise, generate it in each realization

#### starting new section to use setup components
#### THESE NEED TO BE ADJUSTED TO multsame2 - that is, use a new name... the current version has the names from setup

#' @param efthr3 (estinfo) the threshold effect size for communicating about management
#' @param efmn3 (estinfo) the underlying mean change in establishment probability (as a proportion)
#' @param efsem3 (estinfo) the standard deviation of the effect
#' @param extx3 (genlocs) range of x coordinates
#' @param exty3 (genlocs) range of y coordinates
#' @param nodnum3 (genlocs) the number of nodes
#' @param rand3 (genlocs) if TRUE then locations are randomly generated
#' @param loctypei3 info (initvals) the type of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence
#' @param numiniti3 info (initvals) the number of initial locations for presence
#' @param numorpropi3 info (initvals) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param propiniti3 info (initvals) the proportion of initial locations for presence
#' @param loctypee3 estab (initvals) the type of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence
#' @param numinite3 estab (initvals) the number of initial locations for presence
#' @param numorprope3 estab (initvals) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param propinite3 estab (initvals) the proportion of initial locations for presence

#### ending section of startup components that is to be renamed


#' @param distfcn3 com (genmovnet) the function of distance used to estimate movement probability - 'random' (not related to distance) or 'powerlaw' or 'exp' or 'inveuc' (inverse of Euclidean distance)
#' @param expparcn3 com (genmovnet) exponential parameter in XXXXX
#' @param lktypecn3 com (genmovnet) link type, pa is presence/absence (unweighted, occurence/non-occurence), pr is a probability of occurence, wtd1 is general weight
#' @param placn3 com (genmovnet) power law parameter a in ad^(-b)
#' @param plbcn3 com (genmovnet) power law parameter b in ad^(-b)
#' @param tlinkcn3 com (genmovnet) threshold for whether link exists (probability in some cases)

#' @param distfdn3 disp (genmovnet) the function of distance used to estimate movement probability - 'random' (not related to distance) or 'powerlaw' or 'exp' or 'inveuc' (inverse of Euclidean distance)
#' @param exppardn3 disp (genmovnet) exponential parameter in XXXXX
#' @param lktypedn3 disp (genmovnet) link type, pa is presence/absence (unweighted, occurence/non-occurence), pr is a probability of occurence, wtd1 is general weight
#' @param pladn3 disp (genmovnet) power law parameter a in ad^(-b)
#' @param plbdn3 disp (genmovnet) power law parameter b in ad^(-b)
#' @param tlinkdn3 disp (genmovnet) threshold for whether link exists (probability in some cases)

#' @param diag1cn3 com (spreadstep) if T, diagonal of adj mat is given value 1
#' @param mattypecn3 com (spreadstep) 'fromto' indicates that rows are sources and columns are sinks, 'tofrom' indicates that columns are sources and rows are sinks
#' @param max1cn3 com (spreadstep) if T, values in vect2 greater than 1 are replaced by 1
#' @param vect1cn3 com (spreadstep) status of nodes before spread

#' @param diag1dn3 disp (spreadstep) if T, diagonal of adj mat is given value 1
#' @param mattypedn3 disp (spreadstep) 'fromto' indicates that rows are sources and columns are sinks, 'tofrom' indicates that columns are sources and rows are sinks
#' @param max1dn3 disp (spreadstep) if T, values in vect2 greater than 1 are replaced by 1
#' @param vect1dn3 disp (spreadstep) status of nodes before spread

#' @param adpmn3 (makedec) mean probability of adopting management if informed
#' @param adpsdn3 (makedec) sd in truncated normal distribution of probability of adoption

#' @param estpmn3 (estab) mean probability of establishment (new or CONTINUED) in absence of management
#' @param maneffn3 (estab) management effect (proportional reduction in estabp) THIS SHOULD BE GENERATED RATHER THAN SUPPLIED AS INPUT - 2019 removed, should use efmn3, instead
#' @param estpsdn3 (estab) sd in truncated normal distribution

#' @param plotmpn3 if true plots of resulting presence of information and species are generated

#' @keywords simulation experiments
#' @export 
#' @examples
#' x1 <- multdif()
#' TO UPDATE x2 <- multdif(nsims3=2, nsteps3=2, xymatinn=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3),byrow=T,ncol=2), vect1cn=c(1,1,1,0,0,0), vect1dn=c(0,0,0,1,1,1),readxymat3=T)
#' x3 <- multdif(nsims3=2, nsteps3=2, readxymat3=F) # generate everything
#' library(truncnorm)
#' baselinep1 <- multdif(nsims3=20, usethresh3=T, nsteps3=5, readxymat3=F, readcomamn3=F, readdispamn3=F, efmn3=c(0.7), efsem3=c(0.2), efsd3=c(0.2), efthr3=0.5, extx3=c(0,50), exty3=c(0,50), nodenum3=100, rand3=TRUE, loctypei3='random', numorpropi3='prop', propiniti3=0.05, loctypee3='upedge', numorprope3='prop', propinite3=0.05, comip3=c(0.05), spip3=c(0.05), distfcn3='powerlaw', lktypecn3='pr', placn3=c(1), plbcn3=c(0.5), tlinkcn3=0.1, distfdn3='powerlaw', lktypedn3='pr', pladn3=c(1), plbdn3=c(0.5), tlinkdn3=0.1, diag1cn3=T,  max1cn3=T, diag1dn3=T, max1dn3=T, adpmn3=c(0.5), adpsdn3=c(0.2), estpmn3=c(0.5), estpsdn3=c(0.2), plotmpn3=F)
#' baselinep1$multout
#' sens.adpm <- multdif(nsims3=20, usethresh3=T, nsteps3=5, readxymat3=F, readcomamn3=F, readdispamn3=F, efmn3=c(0.7), efsem3=c(0.2), efsd3=c(0.2), efthr3=0.5, extx3=c(0,50), exty3=c(0,50), nodenum3=100, rand3=TRUE, loctypei3='random', numorpropi3='prop', propiniti3=0.05, loctypee3='upedge', numorprope3='prop', propinite3=0.05, comip3=c(0.05), spip3=c(0.05), distfcn3='powerlaw', lktypecn3='pr', placn3=c(1), plbcn3=c(0.5), tlinkcn3=0.1, distfdn3='powerlaw', lktypedn3='pr', pladn3=c(1), plbdn3=c(0.5), tlinkdn3=0.1, diag1cn3=T,  max1cn3=T, diag1dn3=T, max1dn3=T, adpmn3=seq(0,1,0.1), adpsdn3=c(0.2), estpmn3=c(0.5), estpsdn3=c(0.2), plotmpn3=F)
#' jt <- sens.adpm$multout
#' plot(jt$adpmn3, jt$mestab, xlab='Mean probability of adopting management if informed', ylab='Invasive establishment rate (phrasing to be improved)')
#' plot(jt$adpmn3, jt$mdec, xlab='Mean probability of adopting management if informed', ylab='Tech adoption rate (phrasing to be improved)')
#' sens.plbc <- multdif(nsims3=20, usethresh3=T, nsteps3=5, readxymat3=F, readcomamn3=F, readdispamn3=F, efmn3=c(0.7), efsem3=c(0.2), efsd3=c(0.2), efthr3=0.5, extx3=c(0,50), exty3=c(0,50), nodenum3=100, rand3=TRUE, loctypei3='random', numorpropi3='prop', propiniti3=0.05, loctypee3='upedge', numorprope3='prop', propinite3=0.05, comip3=c(0.05), spip3=c(0.05), distfcn3='powerlaw', lktypecn3='pr', placn3=c(1), plbcn3=seq(0,2,0.1), tlinkcn3=0.1, distfdn3='powerlaw', lktypedn3='pr', pladn3=c(1), plbdn3=c(0.5), tlinkdn3=0.1, diag1cn3=T,  max1cn3=T, diag1dn3=T, max1dn3=T, adpmn3=c(0.5), adpsdn3=c(0.2), estpmn3=c(0.5), estpsdn3=c(0.2), plotmpn3=F)
#' jt2 <- sens.plbc$multout
#' plot(jt2$plbcn3, jt2$mestab, xlab='Power law parameter b in ad^(-b)', ylab='Invasive establishment rate (phrasing to be improved)')
#' plot(jt2$plbcn3, jt2$mdec, xlab='Power law parameter b in ad^(-b)', ylab='Tech adoption rate (phrasing to be improved)')
#' plot(jt2$plbcn3, jt2$mcom, xlab='Power law parameter b in ad^(-b)', ylab='Rate of receipt of communication (phrasing to be improved)')


# to do - general test 17
# to do - VERY IMPORTANT for baselinep1 and sensitivity analyses, check the roles of lktypecn3, tlinkcn3, tlinkdn3, diag1cn3, max1cn3, diag1dn3, max1dn3... currently using defaults

# to do - GENERAL TESTING
# to do - confirm that maneffn3 is not doing anything (as it should not be) and remove it
# to do - consider making the type of plot below a function (flexplot) - and an option to provide just one example
# to do - consider how best to use potential outcome of not crossing threshold so that communication occurs...
# to do - add time step number label to figures
# to do - think about in what ways comam and dispam might change from one time step to another...
# to do - IMPORTANT!!! dealing with decision across time... !!!
# to do - IMPORTANT - this has decisions being made again in each time step - not really what was intended, necessarily... need to adjust or decide to keep IMPORTANT - because of nature of onetstep

# to do - note that plots even when plotting options set to F
# to do - there are two variables in conflict - efmn is the mean effect, which should be input.  maneffn3 doesn't make sense as input - it should be generated using efmn
# to do - add mean disp to output
# 2018-07-22 - to do - seems that comip3 and spip3 aren't actually used


multdif <- function(nsims3=30, usethresh3=F, nsteps3=5, xymatinn3, readxymat3=F, vect1cn3=1, vect1dn3=1, comamn3=1, dispamn3=1, readcomamn3=F, readdispamn3=F, 
efmn3=c(0.1,0.8), efsem3=c(0.1), efsd3=c(0.1), efthr3=0.5, 
extx3=c(0,50), exty3=c(0,50), nodenum3=100, rand3=TRUE, loctypei3='random', numiniti3=5, numorpropi3='num', propiniti3=0.05, loctypee3='upedge', numinite3=5, numorprope3='num', propinite3=0.05,
comip3=c(0.05), spip3=c(0.05), distfcn3='powerlaw', expparcn3=1, lktypecn3='pr', placn3=c(1,2), plbcn3=c(1,2), tlinkcn3=0.1, distfdn3='powerlaw', exppardn3=1, lktypedn3='pr', pladn3=c(1,2), plbdn3=c(1,2), tlinkdn3=0.1, diag1cn3=T, mattypecn3='fromto', max1cn3=T, diag1dn3=T, mattypedn3='fromto', max1dn3=T, adpmn3=c(0.1,0.9), adpsdn3=c(0.1), estpmn3=c(0.1,0.9), estpsdn3=c(0.1), plotmpn3=F){

# note that the following line assumes use of powerlaw
# 13 parameters being considered

  ncombs <- length(efmn3)*length(efsem3)*length(efsd3)*length(efthr3)*length(comip3)*length(spip3)*length(placn3)*length(plbcn3)*length(pladn3)*length(plbdn3)*length(adpmn3)*length(adpsdn3)*length(estpmn3)*length(estpsdn3)*length(nodenum3)*length(propiniti3)*length(propinite3)

# ncol depends on how many variables may have multiple values, plus 6 output variables

  vn <- c('efmn3', 'efsem3', 'efsd3', 'efthr3', 'comip3', 'spip3', 'placn3', 'plbcn3', 'pladn3', 'plbdn3', 'adpmn3', 'adpsdn3', 'estpmn3', 'estpsdn3', 'nodenum3', 'propiniti3', 'propinite3', 'mcom', 'mdec', 'mestab', 'vcom', 'vdec', 'vestab')

  multout <- data.frame(matrix(-1,ncol=length(vn),nrow=ncombs))
  names(multout) <- vn
  
  multdetails <- as.list(1:(ncombs*nsims3))

  rowcount <- 1

  for(j1 in 1:length(efmn3)) {
    tefmn3 <- efmn3[j1]
  for(j2 in 1:length(efsem3)) {
    tefsem3 <- efsem3[j2]
  for(j3 in 1:length(efsd3)) {
    tefsd3 <- efsd3[j3]
  for(j3b in 1:length(efthr3)) { # index j3b
    tefthr3 <- efthr3[j3b]
  for(j4 in 1:length(comip3)) {
    tcomip3 <- comip3[j4]
  for(j5 in 1:length(spip3)) {
    tspip3 <- spip3[j5]
  for(j6 in 1:length(placn3)) {
    tplacn3 <- placn3[j6]
  for(j7 in 1:length(plbcn3)) {
    tplbcn3 <- plbcn3[j7]
  for(j8 in 1:length(pladn3)) {
    tpladn3 <- pladn3[j8]
  for(j9 in 1:length(plbdn3)) {
    tplbdn3 <- plbdn3[j9]
  for(j10 in 1:length(adpmn3)) {
    tadpmn3 <- adpmn3[j10]
  for(j11 in 1:length(adpsdn3)) {
    tadpsdn3 <- adpsdn3[j11]
  for(j12 in 1:length(estpmn3)) {
    testpmn3 <- estpmn3[j12]
  for(j13 in 1:length(estpsdn3)) {
    testpsdn3 <- estpsdn3[j13]
  for(j14 in 1:length(nodenum3)) {
    tnodenum3 <- nodenum3[j14]
  for(j15 in 1:length(propiniti3)) {
    tpropiniti3 <- propiniti3[j15]
  for(j16 in 1:length(propinite3)) {
    tpropinite3 <- propinite3[j16]



##### !!!!! to be added in multsame: efmn3, efsem3, efsd3, comip3, spip3, using _tefmn3_ etc in call below

####### NEW VARIABLES NEED THE t PREFIX ADDED WHERE APPROPRIATE
    temp <- multsame2(nsims=nsims3, usethresh=usethresh3, nsteps2=nsteps3, xymatinn2=xymatinn3, readxymat2=readxymat3, vect1cn=vect1cn3, vect1dn=vect1dn3, comamn2=comamn3, dispamn2=dispamn3, readcomamn2=readcomamn3, readdispamn2=readdispamn3, 
efmn22=tefmn3, efsem22=tefsem3, efthr22=tefthr3, extx22=extx3, exty22=exty3, nodenum22=tnodenum3, rand22=rand3, loctypei2=loctypei3, numiniti2=numiniti3, numorpropi2=numorpropi3, propiniti2=tpropiniti3, loctypee2=loctypee3, numinite2=numinite3, numorprope2=numorprope3, propinite2=tpropinite3,
distfcn2=distfcn3, expparcn2=expparcn3, lktypecn2=lktypecn3, placn2=tplacn3, plbcn2=tplbcn3, tlinkcn2=tlinkcn3, distfdn2=distfdn3, exppardn2=exppardn3, lktypedn2=lktypedn3, pladn2=tpladn3, plbdn2=tplbdn3, tlinkdn2=tlinkdn3, diag1cn2=diag1cn3, mattypecn2=mattypecn3, max1cn2=max1cn3, diag1dn2=diag1dn3, mattypedn2=mattypedn3, max1dn2=max1dn3, adpmn2=tadpmn3, adpsdn2=tadpsdn3, estpmn2=testpmn3, estpsdn2=testpsdn3, plotmpn2=plotmpn3)

  multdetails[[rowcount]] <- temp

  multout$efmn3[rowcount] <- tefmn3
  multout$efsem3[rowcount] <- tefsem3
  multout$efsd3[rowcount] <- tefsd3
  multout$efthr3[rowcount] <- tefthr3
  multout$comip3[rowcount] <- tcomip3
  multout$spip3[rowcount] <- tspip3
  multout$placn3[rowcount] <- tplacn3
  multout$plbcn3[rowcount] <- tplbcn3
  multout$pladn3[rowcount] <- tpladn3
  multout$plbdn3[rowcount] <- tplbdn3
  multout$adpmn3[rowcount] <- tadpmn3
  multout$adpsdn3[rowcount] <- tadpsdn3
  multout$estpmn3[rowcount] <- testpmn3
  multout$estpsdn3[rowcount] <- testpsdn3
  multout$nodenum3[rowcount] <- tnodenum3
  multout$propiniti3[rowcount] <- tpropiniti3
  multout$propinite3[rowcount] <- tpropinite3

  multout$mcom[rowcount] <- temp$meancom
  multout$mdec[rowcount] <- temp$meandec
  multout$mestab[rowcount] <- temp$meanestab
  multout$vcom[rowcount] <- var(temp$setocom)
  multout$vdec[rowcount] <- var(temp$setodec)
  multout$vestab[rowcount] <- var(temp$setoestab)

  rowcount <- rowcount + 1
  print(paste("rowcount",rowcount))

  }}}}}}}}}}}}}}}}}

  list(multdetails=multdetails, multout=multout)
}