
#' Estimates the management effect on the probability of establishment of a bioentity, and whether the effect is large enough to trigger communication about the management
#'
#'This function adds a 'science of science' component to INA simulations, allowing consideration of whether managements with lower mean effects and/or higher variability in effect (and/or more limited study, reflected in higher standard error of the mean management effect size) may not have regional effects due to lack of communication.

#'The management effect is a function of the mean effect and the standard error of the mean of the effect.  Another input is the minimum size (as a proportion) of the effect necessary to trigger communication.  The output is the observed change in the probability of establishment, as a proportion, and a logic variable indicating whether the effect is greater than the communication threshold.  (Used by function INAscene)

#' If the communication threshold is not a component of simulations, then communication can occur regardless of management effect size
#'
#' Updated 2020-06-02

#' @param efmn the underlying mean change in establishment probability (as a proportion) by the management being considered
#' @param efsd the standard deviation of the management effect on establishment 
#' @param efthr the threshold management effect size for communicating about the management (if this is set to zero, then communication always occurs - there is no minimum management effect size)
#' @param sampeffort sampling effort, where greater samping effort reduces the error in estimating the management effect
#' @keywords information
#' @export 

#' @examples
#' estinfo(efmn=0.5, efsd=0.5, efthr=0.5, sampeffort=1)
#' estinfo(efmn=0.1, efsd=0.5, efthr=0.5, sampeffort=10)
#' estinfo(efmn=0.5, efsd=0.5, efthr=0.5, sampeffort=10)


estinfo <- function(efmn, efsd, efthr, sampeffort){

  obschange <- rtruncnorm(n=1, a=0, b=1, mean=efmn, sd=efsd/sampeffort)

  com.yes <- obschange > efthr

  list(obschange=obschange, com.yes=com.yes)
}