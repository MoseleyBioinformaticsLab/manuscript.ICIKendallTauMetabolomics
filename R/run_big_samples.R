##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param in_data
##' @return
##' @author rmflight
##' @export
run_big_samples <- function(in_data) {
  future::plan(multicore(workers = 10))
  cor = ici_kendalltau(t(in_data$data))
  cor$type = in_data$type
  cor$n_sample = in_data$n_sample
  future::plan(multicore)
  return(cor)
}
