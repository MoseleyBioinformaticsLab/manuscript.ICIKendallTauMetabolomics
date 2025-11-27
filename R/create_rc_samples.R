##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author rmflight
##' @export
create_lc_samples <- function() {

  base_sample = rlnorm(1000, meanlog = 1, sdlog = 0.5)
  rep_data = add_uniform_noise(2, base_sample, 0.2)
  rep_data

}

left_censor_correlate = function(lc_samples, cut_values = seq(0, 1.5, by = 0.1)){
  censor_cor = purrr::map_df(cut_values, function(in_cut){
    tmp_lc = lc_samples
    tmp_lc[tmp_lc < in_cut] = NA
    n_na = sum(is.na(tmp_lc))
    
    ici_cor = ici_kt(tmp_lc[, 1], tmp_lc[, 2], perspective = "global")
    p_cor = cor.test(tmp_lc[, 1], tmp_lc[, 2], use = "pairwise.complete.obs", method = "pearson")
    k_cor = cor.test(tmp_lc[, 1], tmp_lc[, 2], use = "pairwise.complete.obs", method = "kendall")
    tmp_lc[is.na(tmp_lc)] = 0
    p_cor_0 = cor.test(tmp_lc[, 1], tmp_lc[, 2], method = "pearson")
    k_cor_0 = cor.test(tmp_lc[, 1], tmp_lc[, 2], method = "kendall")
    
    data.frame(cor = c(ici_cor["tau"],
                       p_cor$estimate,
                       k_cor$estimate,
                       p_cor_0$estimate,
                       k_cor_0$estimate),
               p_value = c(ici_cor["pvalue"],
                           p_cor$p.value,
                           k_cor$p.value,
                           p_cor_0$p.value,
                           k_cor_0$p.value),
               which = c("ici",
                         "pearson",
                         "kendall",
                         "pearson_0",
                         "kendall_0"),
               n_na = n_na,
               cutoff = in_cut)
    
  })
  censor_cor
}

lt_left_censor_correlate = function(lc_samples, cut_values = seq(0, 1.5, by = 0.1)){
  censor_cor = purrr::map_df(cut_values, function(in_cut){
    tmp_lc = lc_samples
    tmp_lc[tmp_lc < in_cut] = NA
    n_na = sum(is.na(tmp_lc))
    
    tmp_lc = log(tmp_lc)
    
    ici_cor = ici_kt(tmp_lc[, 1], tmp_lc[, 2], perspective = "global")
    p_cor = cor.test(tmp_lc[, 1], tmp_lc[, 2], use = "pairwise.complete.obs", method = "pearson")
    k_cor = cor.test(tmp_lc[, 1], tmp_lc[, 2], use = "pairwise.complete.obs", method = "kendall")
    tmp_lc[is.na(tmp_lc)] = 0
    p_cor_0 = cor.test(tmp_lc[, 1], tmp_lc[, 2], method = "pearson")
    k_cor_0 = cor.test(tmp_lc[, 1], tmp_lc[, 2], method = "kendall")
    
    data.frame(cor = c(ici_cor["tau"],
                       p_cor$estimate,
                       k_cor$estimate,
                       p_cor_0$estimate,
                       k_cor_0$estimate),
               p_value = c(ici_cor["pvalue"],
                           p_cor$p.value,
                           k_cor$p.value,
                           p_cor_0$p.value,
                           k_cor_0$p.value),
               which = c("ici",
                         "pearson",
                         "kendall",
                         "pearson_0",
                         "kendall_0"),
               n_na = n_na,
               cutoff = in_cut)
    
  })
  censor_cor
}

random_censor_correlate = function(lc_samples, n_na = seq(0, 250, 50), nrep = 100){
  censor_cor = furrr::future_map_dfr(n_na, function(in_na){
    tmp_lc = lc_samples
    n_total = length(lc_samples)
    
    purrr::map_df(seq(1, nrep), function(in_rep){
      na_locs = sample(n_total, in_na)
      tmp_lc[na_locs] = NA
      
      ici_cor = ici_kt(tmp_lc[, 1], tmp_lc[, 2], perspective = "global")
      p_cor = cor.test(tmp_lc[, 1], tmp_lc[, 2], use = "pairwise.complete.obs", method = "pearson")
      k_cor = cor.test(tmp_lc[, 1], tmp_lc[, 2], use = "pairwise.complete.obs", method = "kendall")
      tmp_lc[is.na(tmp_lc)] = 0
      p_cor_0 = cor.test(tmp_lc[, 1], tmp_lc[, 2], method = "pearson")
      k_cor_0 = cor.test(tmp_lc[, 1], tmp_lc[, 2], method = "kendall")
      
      tmp_frame = data.frame(cor = c(ici_cor["tau"],
                                     p_cor$estimate,
                                     k_cor$estimate,
                                     p_cor_0$estimate,
                                     k_cor_0$estimate),
                             p_value = c(ici_cor["pvalue"],
                                         p_cor$p.value,
                                         k_cor$p.value,
                                         p_cor_0$p.value,
                                         k_cor_0$p.value),
                             which = c("ici",
                                       "pearson",
                                       "kendall",
                                       "pearson_0",
                                       "kendall_0"),
                             n_na = in_na,
                             rep = in_rep,
                             cutoff = 0)
      #tmp_frame$na_locs = list(na_locs)
      tmp_frame
    })
  })
  censor_cor
}

add_uniform_noise <- function(n_rep, value, sd, use_zero = FALSE){
  n_value <- length(value)
  
  n_sd <- n_rep * n_value
  
  out_sd <- rnorm(n_sd, 0, sd)
  out_sd <- matrix(out_sd, nrow = n_value, ncol = n_rep)
  
  if (!use_zero){
    tmp_value <- matrix(value, nrow = n_value, ncol = n_rep, byrow = FALSE)
    out_value <- tmp_value + out_sd
  } else {
    out_value <- out_sd
  }
  
  return(out_value)
}
