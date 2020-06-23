
#' Generate geographic locations of nodes in network
#'
#' This function generates the geographic locations of nodes in network, for simulation studies where geographic locations are not observed.
#'
#' Updated 2020-06-02

#' @param extx range of x coordinates
#' @param exty range of y coordinates
#' @param nn the number of nodes
#' @param rand if TRUE then locations are randomly generated (other options for distributing locations will be added)
#' @keywords geography locations
#' @export

#' @examples
#' xymat99 <- genlocs(extx=c(0,50), exty=c(0,50), nn=100, rand=TRUE)


genlocs <- function(extx, exty, nn, rand=TRUE){

  if (rand) {

    xvec <- extx[1] + (extx[2] - extx[1])*runif(n=nn)

    yvec <- exty[1] + (exty[2] - exty[1])*runif(n=nn)

  }

  xyob <- cbind(xvec,yvec)
  xyob
}