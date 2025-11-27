##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param n
##' @return
##' @author rmflight
##' @export
create_sample <- function(n = 1000) {

  #set.seed(1234)
  base_sample = sort(rlnorm(n, meanlog = 1, sdlog = 0.5))
  base_sample

}
