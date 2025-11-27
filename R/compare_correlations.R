##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param where_na
##' @return
##' @author rmflight
##' @export
compare_positive_kt <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  #prog_where = knitrProgressBar::progress_estimated(length(where_na))
  tmp = furrr::future_imap_dfr(where_na, function(use_na, i_na){
    #message(.y)
    #knitrProgressBar::update_progress(prog_where)
    tmp_x = x
    tmp_y = y
    
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    # switching to only introducing NA values in either X or Y,
    # not both, because both makes it NOT symmetric, there is
    # just no way to do it.
    out_val_y = ici_kt(x, tmp_y, perspective = perspective)[[1]]
    out_val_x = ici_kt(tmp_x, y, perspective = perspective)[[1]]
    data.frame(cor = c(out_val_y, out_val_x), i_na = i_na, x_na = c(0, length(x_na)), y_na = c(length(y_na), 0))
  }, .progress = TRUE)
  
  tmp
}

compare_positive_kt_c <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  #prog_where = knitrProgressBar::progress_estimated(length(where_na))
  tmp = furrr::future_map_dbl(where_na, function(use_na){
  #purrr::imap_dbl(where_na, function(use_na, .y){
    #message(.y)
    #knitrProgressBar::update_progress(prog_where)
    tmp_x = x
    tmp_y = y
    
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    ici_kt(tmp_x, tmp_y, perspective = perspective)[[1]]
  #})
  }, .progress = TRUE)
  
  tmp
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param where_na
##' @return
##' @author rmflight
##' @export
compare_negative_kt <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  furrr::future_imap_dfr(where_na, function(use_na, i_na){
    #message(use_na)
    tmp_x = x
    tmp_y = y
    y_na = (n_entry + 1) - (use_na[use_na > n_entry] - n_entry)
    x_na = use_na[use_na <= n_entry]
    
    if ((n_entry == 10)) {
      if (low_indices) {
        y_na = y_na[y_na <= 5]
        y_na = y_na + 5
        x_na = x_na[x_na <= 5]
      } 
    }
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    out_val_y = ici_kt(x, tmp_y, perspective = perspective)[[1]]
    out_val_x = ici_kt(tmp_x, y, perspective = perspective)[[1]]
    data.frame(cor = c(out_val_y, out_val_x), i_na = i_na, x_na = c(0, length(x_na)), y_na = c(length(y_na), 0))
  })
#  , .progress = TRUE)
  
  
}

compare_negative_kt_c <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  furrr::future_map_dbl(where_na, function(use_na){
    #message(.y)
    tmp_x = x
    tmp_y = y
    y_na = use_na[use_na > n_entry] - n_entry + 500
    x_na = use_na[use_na <= n_entry]
    
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    y_na = n_entry - y_na + 1
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    ici_kt(tmp_x, tmp_y, perspective = perspective)[[1]]
  }, .progress = TRUE)
  
  
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param where_na
##' @return
##' @author rmflight
##' @export
compare_negative_pearson <- function(x, y, where_na, low_indices = FALSE, method = "pearson") {
  n_entry = length(x)
  furrr::future_imap_dfr(where_na, function(use_na, i_na){
  #purrr::imap_dfr(where_na, function(use_na, i_na){
    #message(i_na)
    #knitrProgressBar::update_progress(prog_where)
    tmp_x = x
    tmp_y = y
    
    y_na = (n_entry + 1) - (use_na[use_na > n_entry] - n_entry)
    x_na = use_na[use_na <= n_entry]
    
    if ((n_entry == 10)) {
      if (low_indices) {
        y_na = y_na[y_na <= 5]
        y_na = y_na + 5
        x_na = x_na[x_na <= 5]
      } 
    } 

    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    #in_matrix = cbind(tmp_x, tmp_y)
    out_val_y = cor(x, tmp_y, use = "pairwise.complete.obs", method = method)
    out_val_x = cor(tmp_x, y, use = "pairwise.complete.obs", method = method)
    data.frame(cor = c(out_val_y, out_val_x), i_na = i_na, x_na = c(0, length(x_na)), y_na = c(length(y_na), 0))
  })
  #}, .progress = TRUE)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param where_na
##' @return
##' @author rmflight
##' @export
compare_positive_pearson <- function(x, y, where_na, low_indices = FALSE, method = "pearson") {
  n_entry = length(x)
  #prog_where = knitrProgressBar::progress_estimated(length(where_na))
  furrr::future_imap_dfr(where_na, function(use_na, i_na){
    #message(.y)
    #knitrProgressBar::update_progress(prog_where)
    tmp_x = x
    tmp_y = y
    
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    in_matrix = cbind(tmp_x, tmp_y)
    out_val_y = cor(x, tmp_y, use = "pairwise.complete.obs", method = method)
    out_val_x = cor(tmp_x, y, use = "pairwise.complete.obs", method = method)
    data.frame(cor = c(out_val_y, out_val_x), i_na = i_na, x_na = c(0, length(x_na)), y_na = c(length(y_na), 0))
  }, .progress = TRUE)
  
  
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param where_na
##' @return
##' @author rmflight
##' @export
compare_negative_kendall <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  furrr::future_imap_dfr(where_na, function(use_na, i_na){
    #purrr::imap_dfr(where_na, function(use_na, i_na){
    #message(i_na)
    #knitrProgressBar::update_progress(prog_where)
    tmp_x = x
    tmp_y = y
    
    y_na = (n_entry + 1) - (use_na[use_na > n_entry] - n_entry)
    x_na = use_na[use_na <= n_entry]
    
    if ((n_entry == 10)) {
      if (low_indices) {
        y_na = y_na[y_na <= 5]
        y_na = y_na + 5
        x_na = x_na[x_na <= 5]
      } 
    } 
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    #in_matrix = cbind(tmp_x, tmp_y)
    
    
    out_val_y = ici_kt(x[-y_na], tmp_y[-y_na], perspective = perspective)[1]
    out_val_x = ici_kt(tmp_x[-x_na], y[-x_na], perspective = perspective)[1]
    data.frame(cor = c(out_val_y, out_val_x), i_na = i_na, x_na = c(0, length(x_na)), y_na = c(length(y_na), 0))
  })
  #}, .progress = TRUE)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param where_na
##' @return
##' @author rmflight
##' @export
compare_positive_kendall <- function(x, y, where_na, low_indices = FALSE, perspective = "global") {
  n_entry = length(x)
  #prog_where = knitrProgressBar::progress_estimated(length(where_na))
  furrr::future_imap_dfr(where_na, function(use_na, i_na){
    #message(.y)
    #knitrProgressBar::update_progress(prog_where)
    tmp_x = x
    tmp_y = y
    
    y_na = use_na[use_na > n_entry] - n_entry
    x_na = use_na[use_na <= n_entry]
    
    if (low_indices) {
      y_na = y_na[y_na <= 5]
      x_na = x_na[x_na <= 5]
    }
    
    
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    in_matrix = cbind(tmp_x, tmp_y)
    out_val_y = ici_kt(x[-y_na], tmp_y[-y_na], perspective = perspective)[1]
    out_val_x = ici_kt(tmp_x[-x_na], y[-x_na], perspective = perspective)[1]
    data.frame(cor = c(out_val_y, out_val_x), i_na = i_na, x_na = c(0, length(x_na)), y_na = c(length(y_na), 0))
  }, .progress = TRUE)
  
  
}
