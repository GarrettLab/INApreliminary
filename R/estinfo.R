
#' Estimates the management effect on the probability of establishment of a bioentity, and whether the effect is large enough to trigger communication about the management
#'
#'This function adds a 'science of science' component to INA simulations, allowing consideration of whether managements with lower mean effects and/or higher variability in effect (and/or more limited study, reflected in higher standard error of the mean management effect size) may not have regional effects due to lack of communication.

#'The management effect is a function of the mean effect and the standard error of the mean of the effect.  Another input is the minimum size (as a proportion) of the effect necessary to trigger communication.  The output is the observed change in the probability of establishment, as a proportion, and a logic variable indicating whether the effect is greater than the communication threshold.  (Used by function INAscene, called directly by function setup2.)

#' If the communication threshold is not a component of simulations, then communication can occur regardless of management effect size
#'
#' Updated 2020-09-05

#' @param maneffmean4s the underlying mean change in establishment probability (as a proportion) by the management being considered
#' @param maneffsd4s the standard deviation of the management effect on establishment 
#' @param maneffthresh4 the threshold management effect size for communicating about the management (if this is set to zero, then communication always occurs - there is no minimum management effect size)
#' @param sampeffort4 sampling effort, where greater samping effort reduces the error in estimating the management effect
#' @keywords information
#' @export 
#' @import truncnorm

#' @examples

#' estinfo(maneffmean4s=0.5, maneffsd4s=0.5, maneffthresh4=0.5, sampeffort4=1)

#' estinfo(maneffmean4s=0.1, maneffsd4s=0.5, maneffthresh4=0.5, sampeffort4=10)

#' estinfo(maneffmean4s=0.5, maneffsd4s=0.5, maneffthresh4=0.5, sampeffort4=10)


estinfo <- function(maneffmean4s, maneffsd4s, maneffthresh4, sampeffort4){

  obschange <- sum(rtruncnorm(n=sampeffort4, a=0, b=1, mean=maneffmean4s, sd=maneffsd4s))/sampeffort4

  com.yes <- obschange > maneffthresh4

  list(obschange=obschange, com.yes=com.yes)
}