##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param x
##' @return
##' @author rmflight
##' @export
create_na_indices <- function(x = 20) {

  n_na = seq(1, x)
  
  where_na = purrr::map(n_na, function(in_na){
    na_comb = combn(x, in_na)
    asplit(na_comb, 2)
  })
  where_na = unlist(where_na, recursive = FALSE)

}

create_random_na = function(x = 500, n_by = 10, nrep = 5, xloc = seq(1, 500), yloc = seq(1001, 1500)) {
  n_na = seq(1, x, by = n_by)
  where_na = vector("list", length = (length(n_na)^2 * nrep^2))
  save_loc = 1
  for (x_na in n_na) {
    for (xrep in seq(1, nrep)) {
      for (y_na in n_na) {
        for (yrep in seq(1, nrep)) {
          where_na[[save_loc]] = c(sample(xloc, x_na), sample(yloc, y_na))
          save_loc = save_loc + 1
          #message(save_loc)
        }
      }
    }
  }
  where_na
}
