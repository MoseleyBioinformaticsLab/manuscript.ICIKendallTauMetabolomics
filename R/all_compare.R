##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param where_na
##' @return
##' @author rmflight
##' @export
all_kendalltau <- function(s1, s2, sneg, where_na, perspective = "global") {
  n_entry = length(s1)
  #prog_where = knitrProgressBar::progress_estimated(length(where_na))
  tmp = furrr::future_map_dfr(where_na, function(use_na){
    #message(.y)
    #knitrProgressBar::update_progress(prog_where)
    tmp_s1 = s1
    tmp_s2 = s2
    tmp_neg = sneg
    
    s2_na = use_na[use_na > n_entry] - n_entry
    s1_na = use_na[use_na <= n_entry]
    neg_na = (n_entry + 1) - (use_na[use_na > n_entry] - n_entry)
    
    tmp_s1[s1_na] = NA
    tmp_s2[s2_na] = NA
    tmp_neg[neg_na] = NA
    s1_na = sum(is.na(tmp_s1))
    s2_na = sum(is.na(tmp_s2))
    neg_na = sum(is.na(tmp_neg))
    
    s1_s2 = sum(is.na(tmp_s1) & is.na(tmp_s2))
    s1_neg = sum(is.na(tmp_s1) & is.na(tmp_neg))
    ici_s1 = ici_kt(tmp_s1, tmp_s2, perspective = perspective)[[1]]
    ici_neg = ici_kt(tmp_s1, tmp_neg, perspective = perspective)[[1]]
    
    kt_s1 = tmp_s1
    kt_s1[is.na(kt_s1)] = 0
    kt_s2 = tmp_s2
    kt_s2[is.na(kt_s2)] = 0
    kt_neg = tmp_neg
    kt_neg[is.na(kt_neg)] = 0
    
    ktcor_s1 = cor(kt_s1, kt_s2, method = "kendall")
    ktcor_neg = cor(kt_s1, kt_neg, method = "kendall")
    
    data.frame(ici_kt = c(ici_s1, ici_neg),
               kendall = c(ktcor_s1, ktcor_neg),
               x_na = c(s1_na, s1_na),
               y_na = c(s2_na, neg_na),
               c_na = c(s1_s2, s1_neg),
               comp = c("Positive", "Negative"))
  })
  tmp
}

