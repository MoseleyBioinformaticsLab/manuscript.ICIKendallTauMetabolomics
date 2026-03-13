##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param select_ss_small
##' @return
##' @author rmflight
##' @export
run_small_samples <- function(in_data) {

  future::plan(multicore(workers = 1))
  cor = ici_kendalltau(t(in_data$data))
  cor$type = in_data$type
  cor$n_sample = in_data$n_sample
  cor
  future::plan(multicore)
  
  return(cor)

}
