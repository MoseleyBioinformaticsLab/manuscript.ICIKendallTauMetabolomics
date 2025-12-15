##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param n
##' @return
##' @author rmflight
##' @export
create_normal_sample <- function(n = 1000, seed) {
  #set.seed(1234)
  base_sample = withr::with_seed(seed = seed, code = {
    sort(rlnorm(n, meanlog = 1, sdlog = 0.5))
  })

  base_sample
}

create_unif_noise_sample = function(base_sample, seed) {
  n = length(base_sample)
  noise_sample = withr::with_seed(seed = seed, code = {
    runif(n, -0.5, 0.5)
  })

  out_sample = (base_sample + noise_sample) |> sort()

  out_sample
}

create_outlier_sample = function(base_sample, perc = 0.005, seed) {
  n1 = length(base_sample)
  n2 = length(base_sample) * perc

  noise_sample = withr::with_seed(seed = seed, code = {
    unif_noise = runif(n = n1, -0.5, 0.5)
    log_noise = sort(rlnorm(n = n2, meanlog = 1.2, sdlog = 0.1))
    unif_noise[seq((n1 - n2) + 1, n1)] = unif_noise[seq((n1 - n2) + 1, n1)] +
      log_noise
    unif_noise
  })

  out_sample = (base_sample + noise_sample) |> sort()
  out_sample
}
