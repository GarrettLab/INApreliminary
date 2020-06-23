
#' Evaluates scenarios in an impact network analysis (INA)
#'
#' This function implements and summarizes multiple simulations across a designated range of parameter values
#'
#' Updated 2020-06-21

#' @param nreals number of realizations to be evaluated
#' @param ntsteps number of time steps to be evaluated
#' @param doplot if true, plots of resulting presence of information and bioentity are generated

#' @param readgeocoords if T, read in geocoords - otherwise, generate it in each realization
#' @param readorgengeo if 'read', read in geocoords - if 'gen', generate it in each realization (redundant wtih readgeocoords)
#' @param geocoords matrix of x,y coordinates of nodes

#' @param numnodes (genlocs:nodenum) the number of nodes, can be a vector of different numbers of nodes for scenario comparisons
#' @param xrange (genlocs:extx) range of x coordinates, e.g., c(0,50)
#' @param yrange (genlocs:exty) range of y coordinates, e.g., c(0,20)
#' @param randgeo (genlocs:rand) if TRUE then locations are randomly generated (the only location generation option for the moment)

#' @param readinitinfo info (setup2) if T, the initial values for the vector of starting locations for the presence of information are read in rather than generated
#' @param initinfo info (setup2) the vector of initial values read in if readinitinfo == T, a vector with lenth equal to the number of nodes and entries 1s or 0s with 1s indicating the initial presence of information

#' @param initinfo.norp info (initvals:numorpropi) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param initinfo.n info (initvals:numiniti) the number of initial locations for presence
#' @param initinfo.p info (initvals:propiniti) the proportion of initial locations for presence
#' @param initinfo.dist info (initvals:loctypei) the type of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence

#' @param readinitbio bio (setup2) if T, the initial values for the vector of starting locations for the presence of the bioentity are read in rather than generated
#' @param initbio bio (setup2) the vector of initial values read in if readinitbio.s == T, a vector with lenth equal to the number of nodes and entries 1s or 0s with 1s indicating the initial presence of the bioentity

#' @param initbio.norp estab (initvals:numorprope/b) 'num' indicates initial number for presence, 'prop' indicates initial proportion
#' @param initbio.n estab (initvals:numinite/b) the number of initial locations for presence
#' @param initbio.p estab (initvals:propinite/b) the proportion of initial locations for presence
#' @param initbio.dist estab (initvals:loctypee/b) the type of locations where initial presence occurs: 'random' indicates all equally likely, 'upedge' indicates that nodes closest to the upper edge have presence, 'rightedge' indicates that nodes closest to the right edge have presence

#' @param readseam if T, the socioeconomic network adjacency matrix is read in rather than being generated from the outset
#' @param seam socioeconnomic network adjacency matrix, read in if readseam=T

#' @param seamdist com (genmovnet:distfcn) the function of distance used to estimate movement probability - 'random' (not related to distance) or 'powerlaw'
#' @param seamrandp (genmovnet) probability of link existence in random network generated for the socioeconomic network 
#' @param seampla com (genmovnet:placn) power law parameter a in ad^(-b)
#' @param seamplb com (genmovnet:pllbcn) power law parameter b in ad^(-b)
#' @param seamlinktype com (genmovnet:lktypecn) link type, pa is presence/absence (unweighted, occurence/non-occurence), pr is a probability of occurence, wtd1 is general weight
#' @param seamthresh com (genmovnet:tlinkcn) threshold for whether link exists (probability in some cases)

#' @param seamdiag1 com (spreadstep:diag1cn) if T, diagonal of adj mat is given value 1
#' @param seamformat com (spreadstep:mattypecn) 'fromto' indicates that rows are sources and columns are sinks, 'tofrom' indicates that columns are sources and rows are sinks
#' @param seammax1 com (spreadstep:max1cn) if T, values in vect2 greater than 1 are replaced by 1

#' @param readbpam if T, the biophysical network adjacency matrix, describing dispersal likelihoods, is read in rather than being generated from the outset
#' @param bpam biophysical adjacency matrix, read in if readbpam=T

#' @param bpamdist disp (genmovnet:distfdn) the function of distance used to estimate movement probability - 'random' (not related to distance) or 'powerlaw' 
#' @param bpamrandp (genmovnet) if bpamdist='random', the probability a link exists in the biophysical network adjacency matrix
#' @param bpampla disp (genmovnet:pladn) power law parameter a in ad^(-b)
#' @param bpamplb disp (genmovnet:plbdn) power law parameter b in ad^(-b)
#' @param bpamlinktype disp (genmovnet:lktypedn) link type, pa is presence/absence (unweighted, occurence/non-occurence), pr is a probability of occurence, wtd1 is general weight
#' @param bpamthresh disp (genmovnet:tlinkdn) threshold for whether link exists (probability in some cases)

#' @param bpamdiag1 disp (spreadstep:diag1dn) if T, diagonal of adj mat is given value 1
#' @param bpamformat disp (spreadstep:mattypedn) 'fromto' indicates that rows are sources and columns are sinks, 'tofrom' indicates that columns are sources and rows are sinks
#' @param bpammax1 disp (spreadstep:max1dn) if T, values in vect2 greater than 1 are replaced by 1

#' @param readprobadoptvec if T, read in a vector of probabilities of adoption for each node in socioeconomic network - if F, generate this vector
#' @param probadoptvec vector of probabilities of adoption for each node in socioeconomic network, read in if readprobadoptvec=T

#' @param probadoptmean (makedec:adpmn) mean probability of adopting management if informed
#' @param probadoptsd (makedec:adpsdn) sd in truncated normal distribution of probability of adoption

#' @param readprobestabvec if T, read in a vector of probabilities of establishment for each node in biophysical network - if F, generate this vector
#' @param probestabvec vector of probabilities of establishment for each node in biophysical network, read in if readprobestabvec=T

#' @param probestabmean (estab:estpmn) mean probability of establishment (new or CONTINUED) in absence of management
#' @param probestabsd (estab:estpsdn) sd of probability of establishment in absence of management in truncated normal distribution

#' @param maneffmean (estinfo:efmn) the underlying mean change in establishment probability (as a proportion) resulting from management technology
#' @param maneffsd (estinfo) the standard deviation of the management technology effect

#' @param usethreshman if T, the threshold for management is used, so that information is never present anywhere unless the management effect estimate exceeds the threshold
#' @param maneffthresh (estinfo:efthr) the threshold effect size for communicating about management
#' @param sampeffort (estinfo) the sampling effort, where increasing sampling effort results in a more precise estimate of the management effect size




#' @keywords simulation experiments
#' @export 

#' @examples

#' library(truncnorm)


#' j25.readgeocoords <- INAscene(nreals=3, ntsteps=3, doplot=F, readorgengeo='read', readgeocoords=T, geocoords=matrix(c(1,1, 1,2, 1,3, 2,1, 2,2, 2,3),byrow=T,ncol=2), numnodes=NA, xrange=NA, yrange=NA, randgeo=F, readinitinfo=T, initinfo=c(1,1,1,0,0,0), initinfo.norp=NA, initinfo.n=NA, initinfo.p=NA, initinfo.dist=NA, readinitbio=T, initbio=c(0,0,0,1,1,1), initbio.norp=NA, initbio.n=NA, initbio.p=NA,  initbio.dist=NA, readseam=F, seam=NA, seamdist='random', seamrandp=c(0.01,0.05,0.1,0.5), seampla=NA, seamplb=NA, seamlinktype='pa', seamthresh=NA, seamdiag1=T, seamformat='fromto', seammax1=T, readbpam=F, bpam=NA, bpamdist='random', bpamrandp=0.1, bpampla=NA, bpamplb=NA, bpamlinktype='pa', bpamthresh=NA, bpamdiag1=T, bpamformat='fromto', bpammax1=T, readprobadoptvec=F, probadoptvec=NA, probadoptmean=0.1, probadoptsd=0.1, readprobestabvec=F, probestabvec=NA, probestabmean=0.1, probestabsd=0.1, maneffmean=0.5, maneffsd=0.1, usethreshman=F, maneffthresh=NA, sampeffort=NA) 

#' j25.baseline <- INAscene(nreals=3, ntsteps=3, doplot=F, readorgengeo='gen', readgeocoords=F, geocoords=NA, numnodes=50, xrange=c(0,10), yrange=c(0,20), randgeo=T, readinitinfo=F, initinfo=NA, initinfo.norp='num', initinfo.n=5, initinfo.p=NA, initinfo.dist='random', readinitbio=F, initbio=NA, initbio.norp='num', initbio.n=5, initbio.p=NA,  initbio.dist='random', readseam=F, seam=NA, seamdist='powerlaw', seamrandp=NA, seampla=1, seamplb=0.5, seamlinktype='pa', seamthresh=0.1, seamdiag1=T, seamformat='fromto', seammax1=T, readbpam=F, bpam=NA, bpamdist='powerlaw', bpamrandp=NA, bpampla=1, bpamplb=0.5, bpamlinktype='pa', bpamthresh=0.5, bpamdiag1=T, bpamformat='fromto', bpammax1=T, readprobadoptvec=F, probadoptvec=NA, probadoptmean=0.5, probadoptsd=0.2, readprobestabvec=F, probestabvec=NA, probestabmean=0.5, probestabsd=0.2, maneffmean=0.5, maneffsd=0.2, usethreshman=T, maneffthresh=0.5, sampeffort=2) 

#' j25.baseline$multout

#' sens.probadoptmean <- INAscene(nreals=15, ntsteps=3, doplot=F, readorgengeo='gen', readgeocoords=F, geocoords=NA, numnodes=50, xrange=c(0,10), yrange=c(0,20), randgeo=T, readinitinfo=F, initinfo=NA, initinfo.norp='num', initinfo.n=5, initinfo.p=NA, initinfo.dist='random', readinitbio=F, initbio=NA, initbio.norp='num', initbio.n=5, initbio.p=NA,  initbio.dist='random', readseam=F, seam=NA, seamdist='powerlaw', seamrandp=NA, seampla=1, seamplb=0.5, seamlinktype='pa', seamthresh=0.1, seamdiag1=T, seamformat='fromto', seammax1=T, readbpam=F, bpam=NA, bpamdist='powerlaw', bpamrandp=NA, bpampla=1, bpamplb=0.5, bpamlinktype='pa', bpamthresh=0.5, bpamdiag1=T, bpamformat='fromto', bpammax1=T, readprobadoptvec=F, probadoptvec=NA, probadoptmean=seq(0,1,0.1), probadoptsd=0.2, readprobestabvec=F, probestabvec=NA, probestabmean=0.5, probestabsd=0.2, maneffmean=0.9, maneffsd=0.2, usethreshman=T, maneffthresh=0.3, sampeffort=2) 

#' jt <- sens.probadoptmean$multout
#' plot(jt$probadoptmean, jt$mestab, xlab='Mean probability of adopting technology if informed', ylab='Proportion nodes with bioentity')
#' plot(jt$probadoptmean, jt$mdec, xlab='Mean probability of adopting technology if informed', ylab='Proportion nodes with technology adoption')

#' sens.seamplb <- INAscene(nreals=15, ntsteps=3, doplot=F, readorgengeo='gen', readgeocoords=F, geocoords=NA, numnodes=50, xrange=c(0,50), yrange=c(0,50), randgeo=T, readinitinfo=F, initinfo=NA, initinfo.norp='num', initinfo.n=5, initinfo.p=NA, initinfo.dist='random', readinitbio=F, initbio=NA, initbio.norp='num', initbio.n=5, initbio.p=NA,  initbio.dist='random', readseam=F, seam=NA, seamdist='powerlaw', seamrandp=NA, seampla=1, seamplb=seq(0,2,0.1), seamlinktype='pa', seamthresh=0.1, seamdiag1=T, seamformat='fromto', seammax1=T, readbpam=F, bpam=NA, bpamdist='powerlaw', bpamrandp=NA, bpampla=1, bpamplb=0.5, bpamlinktype='pa', bpamthresh=0.5, bpamdiag1=T, bpamformat='fromto', bpammax1=T, readprobadoptvec=F, probadoptvec=NA, probadoptmean=0.7, probadoptsd=0.2, readprobestabvec=F, probestabvec=NA, probestabmean=0.5, probestabsd=0.2, maneffmean=0.9, maneffsd=0.2, usethreshman=T, maneffthresh=0.3, sampeffort=2) 

#' jt2 <- sens.seamplb$multout
#' plot(jt2$seamplb, jt2$mestab, xlab='Power law parameter b in ad^(-b) for communication links', ylab='Proportion nodes with bioentity')
#' plot(jt2$seamplb, jt2$mdec, xlab='Power law parameter b in ad^(-b) for communication links', ylab='Proportion nodes with technology adoption')
#' plot(jt2$seamplb, jt2$mcom, xlab='Power law parameter b in ad^(-b) for communication links', ylab='Proportion nodes with information about technology')


#' j25.seamrandp <- INAscene(nreals=3, ntsteps=3, doplot=F, readorgengeo='gen', readgeocoords=F, geocoords=NA, numnodes=50, xrange=c(0,10), yrange=c(0,20), randgeo=T, readinitinfo=F, initinfo=NA, initinfo.norp='num', initinfo.n=5, initinfo.p=NA, initinfo.dist='random', readinitbio=F, initbio=NA, initbio.norp='num', initbio.n=5, initbio.p=NA,  initbio.dist='random', readseam=F, seam=NA, seamdist='random', seamrandp=c(0.01,0.05,0.1,0.5), seampla=NA, seamplb=NA, seamlinktype='pa', seamthresh=NA, seamdiag1=T, seamformat='fromto', seammax1=T, readbpam=F, bpam=NA, bpamdist='random', bpamrandp=0.1, bpampla=NA, bpamplb=NA, bpamlinktype='pa', bpamthresh=NA, bpamdiag1=T, bpamformat='fromto', bpammax1=T, readprobadoptvec=F, probadoptvec=NA, probadoptmean=0.1, probadoptsd=0.1, readprobestabvec=F, probestabvec=NA, probestabmean=0.1, probestabsd=0.1, maneffmean=0.5, maneffsd=0.1, usethreshman=F, maneffthresh=NA, sampeffort=NA) 




INAscene <- function(nreals, ntsteps, doplot=F, readorgengeo, readgeocoords, geocoords=NA, numnodes=NA, xrange=NA, yrange=NA, randgeo=NA, readinitinfo, initinfo=NA, initinfo.norp=NA, initinfo.n=NA, initinfo.p=NA, initinfo.dist=NA, readinitbio, initbio=NA, initbio.norp=NA, initbio.n=NA, initbio.p=NA,  initbio.dist=NA, readseam, seam=NA, seamdist=NA, seamrandp=NA, seampla=NA, seamplb=NA, seamlinktype=NA, seamthresh=NA, seamdiag1=NA, seamformat=NA, seammax1=NA, readbpam, bpam=NA, bpamdist=NA, bpamrandp=NA, bpampla=NA, bpamplb=NA, bpamlinktype=NA, bpamthresh=NA, bpamdiag1=NA, bpamformat=NA, bpammax1=NA, readprobadoptvec, probadoptvec=NA, probadoptmean=NA, probadoptsd=NA, readprobestabvec, probestabvec=NA, probestabmean=NA, probestabsd=NA, maneffmean=NA, maneffsd=NA, usethreshman, maneffthresh=NA, sampeffort=NA) {

### prepare the output matrix, taking into account which
###   variables are used, as a function of which options
###   were selected

# create initial output matrix, starting with only the 6
#   output variables

vn <- c('mcom', 'mdec', 'mestab', 'vcom', 'vdec', 'vestab')

multout <- data.frame(matrix(-99,ncol=length(vn),nrow=1))

names(multout) <- vn


# function expandoutmat is used to expand the output matrix 
#   (add rows) as more parameter levels are considered

expandoutmat <- function(startoutmat, newvar){
  endmat <- startoutmat
  if (length(newvar) > 1){
    for (i2 in 2:length(newvar)) {
      endmat <- rbind(endmat,startoutmat)
    }
  } 
  endmat
}

# 1. Add to output matrix, taking into account whether
#   reading in geographic locations -or- generating them
# numnodes can be a vector of values

if (readgeocoords==F) {
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(startoutmat=multout,newvar=numnodes)
  multout$numnodes <- rep(numnodes, each=nrowsbefore)
}

# 2.1. Add to output matrix, taking into account whether
#   reading in which nodes have info initially -or- generating
#   the node locations
# initinfo.n or initinfo.p can be a vector of values

if (readinitinfo==F) {
  if (initinfo.norp=='num') {
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=initinfo.n)
    multout$initinfo.n <- rep(initinfo.n, each=nrowsbefore)
  } else if (initinfo.norp=='prop'){
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=initinfo.p)
    multout$initinfo.p <- rep(initinfo.p, each=nrowsbefore)
  }
}

# 2.2. Add to output matrix, taking into account whether
#   reading in which nodes have the bioentity initially -or-
#   generating the node locations
# initbio.n or initbio.p can be a vector of values

if (readinitbio==F) {
  if (initbio.norp=='num') {
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(startoutmat=multout,newvar=initbio.n)
    multout$initbio.n <- rep(initbio.n, each=nrowsbefore)
  } else if (initbio.norp=='prop'){
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=initbio.p)
    multout$initbio.p <- rep(initbio.p, each=nrowsbefore)
  }
}

# 3.1. Add to output matrix, taking into account whether
#   reading in the adjacency matrix for the socioeconomic
#     network -or- generating this adjacency matrix
# seamrandp or seampla and seamplb and seamthresh 
#   can be vectors of values

if (readseam==F) {
  if (seamdist=='random') {
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(startoutmat=multout,newvar=seamrandp)
    multout$seamrandp <- rep(seamrandp, each=nrowsbefore)
  } else if (seamdist=='powerlaw'){
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=seampla)
    multout$seampla <- rep(seampla, each=nrowsbefore)
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=seamplb)
    multout$seamplb <- rep(seamplb, each=nrowsbefore)
  }
  if (seamlinktype=='pa') {
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=seamthresh)
    multout$seamthresh <- rep(seamthresh, each=nrowsbefore)
  }
}

# 3.2. Add to output matrix, taking into account whether
#   reading in the adjacency matrix for the biophysical
#     network -or- generating this adjacency matrix
# bpamrandp or bpampla and bpamplb and bpamthresh 
#   can be vectors of values

if (readbpam==F) {
  if (bpamdist=='random') {
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(startoutmat=multout,newvar=bpamrandp)
    multout$bpamrandp <- rep(bpamrandp, each=nrowsbefore)
  } else if (bpamdist=='powerlaw'){
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=bpampla)
    multout$bpampla <- rep(bpampla, each=nrowsbefore)
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=bpamplb)
    multout$bpamplb <- rep(bpamplb, each=nrowsbefore)
  }
  if (bpamlinktype=='pa') {
    nrowsbefore <- nrow(multout)
    multout <- expandoutmat(multout,newvar=bpamthresh)
    multout$bpamthresh <- rep(bpamthresh, each=nrowsbefore)
  }
}

# 4.1. Add to output matrix, taking into account whether
#   reading in a vector of the probability of adoption
#   -or- generating this vector
# probadoptmean and probadoptsd can be vectors of values

if (readprobadoptvec==F) {
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(multout,newvar=probadoptmean)
  multout$probadoptmean <- rep(probadoptmean, each=nrowsbefore)
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(startoutmat=multout,newvar=probadoptsd)
  multout$probadoptsd <- rep(probadoptsd, each=nrowsbefore)
}

# 4.2. Add to output matrix, taking into account whether
#   reading in a vector of the probability of establishment in
#   the absence of management -or- generating this vector
# probestabmean, probestabsd, maneffmean, maneffsd can be
#   vector of values

if (readprobestabvec==F) {
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(multout,newvar=probestabmean)
  multout$probestabmean <- rep(probestabmean, each=nrowsbefore)
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(startoutmat=multout,newvar=probestabsd)
  multout$probestabsd <- rep(probestabsd, each=nrowsbefore)
}
nrowsbefore <- nrow(multout)
multout <- expandoutmat(startoutmat=multout,newvar=maneffmean)
multout$maneffmean <- rep(maneffmean, each=nrowsbefore)
nrowsbefore <- nrow(multout)
multout <- expandoutmat(startoutmat=multout,newvar=maneffsd)
multout$maneffsd <- rep(maneffsd, each=nrowsbefore)


# 5. Add to output matrix, taking into account whether
#   the estimated management effect must exceed a 
#   threshold for communication to take place
# maneffthresh and sampeffort can be a vector of values

if (usethreshman==TRUE) {
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(multout,newvar=maneffthresh)
  multout$maneffthresh <- rep(maneffthresh, each=nrowsbefore)
  nrowsbefore <- nrow(multout)
  multout <- expandoutmat(multout,newvar=sampeffort)
  multout$sampeffort <- rep(sampeffort, each=nrowsbefore)
}


# number of rows in multout = number of parameter combinations
ncombs <- nrow(multout) 

# number of columns in multout = number of parameters varying 
#   plus 6 response variables
ncmultout <- ncol(multout) 

# set up an object to keep all the output
multdetails <- as.list(1:(ncombs*nreals))

# keep track of which row is being evaluated
rowcount <- 1

### Evaluate the scenario for each parameter combination

for (j22 in 1:ncombs) {

  # prepare input values for multsame2 taking the values
  #   from the current row of multout (params start at 7)
  for (i22 in 7:ncmultout) {
    assign(names(multout)[i22], value=multout[j22,i22])
  }


  temp <- multsame2(nreals=nreals, nsteps2=ntsteps, plotmpn2=doplot, readorgenxy22=readorgengeo, readxymat2=readgeocoords, xymatinn2=geocoords, nodenum22=numnodes, extx22=xrange, exty22=yrange, rand22=randgeo, readinitinfo.s2=readinitinfo, initinfo.s2=initinfo, numorpropi2=initinfo.norp, numiniti2=initinfo.n, propiniti2=initinfo.p, loctypei2=initinfo.dist,  readinitbio.s2=readinitbio,  initbio.s2=initbio,  numorpropb2=initbio.norp, numinitb2=initbio.n, propinitb2=initbio.p,  loctypeb2=initbio.dist, readcomamn2=readseam, comamn2=seam,  distfcn2=seamdist, randpcn2=seamrandp, placn2=seampla, plbcn2=seamplb, lktypecn2=seamlinktype, tlinkcn2=seamthresh, diag1cn2=seamdiag1, mattypecn2=seamformat, max1cn2=seammax1, readdispamn2=readbpam,  dispamn2=bpam, distfdn2=bpamdist, randpdn2=bpamrandp, pladn2=bpampla, plbdn2=bpamplb, lktypedn2=bpamlinktype, tlinkdn2=bpamthresh, diag1dn2=bpamdiag1, mattypedn2=bpamformat, max1dn2=bpammax1, readpadopt2=readprobadoptvec, padopt2=probadoptvec, adpmn2=probadoptmean, adpsdn2=probadoptsd, readpestabn2=readprobestabvec, pestabn2=probestabvec, estpmn2=probestabmean, estpsdn2=probestabsd, efmn22=maneffmean, efsd22=maneffsd, usethresh=usethreshman, efthr22=maneffthresh, sampeffort22=sampeffort)


  multdetails[[rowcount]] <- temp

  multout$mcom[rowcount] <- temp$meancom
  multout$mdec[rowcount] <- temp$meandec
  multout$mestab[rowcount] <- temp$meanestab
  multout$vcom[rowcount] <- var(temp$setocom)
  multout$vdec[rowcount] <- var(temp$setodec)
  multout$vestab[rowcount] <- var(temp$setoestab)

  print(paste("ending row",rowcount))
  rowcount <- rowcount + 1

}

  list(multdetails=multdetails, multout=multout)
}