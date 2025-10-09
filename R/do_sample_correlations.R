median_normalize = function(in_counts) {
  median_counts = apply(in_counts, 2, median, na.rm = TRUE)
  median_counts[median_counts == 0] = 1
  median_matrix = matrix(
    median_counts,
    nrow = nrow(in_counts),
    ncol = ncol(in_counts),
    byrow = TRUE
  )

  norm_counts = in_counts / median_matrix
  return(norm_counts)
}

run_cor_everyway_new = function(id, smd_file, type, value_check) {
  # id = samples_methods$id[1]
  # smd_file = samples_methods$smd_file[1]

  # id = correlation_samples$id[1]
  # smd_file = correlation_samples$smd_file[1]

  sample_smd = readRDS(smd_file)

  sample_counts = assays(sample_smd)$counts
  sample_info = colData(sample_smd) |> as.data.frame()

  rank_data = ICIKendallTau::rank_order_data(
    sample_counts,
    sample_classes = sample_info$factors
  )

  rank_info = purrr::map(rank_data, \(x) {
    x$n_na_rank
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(factors = split)

  censorship_test = ICIKendallTau::test_left_censorship(
    sample_counts,
    sample_classes = sample_info$factors
  )

  keep_counts = keep_non_missing_percentage(
    sample_counts,
    sample_classes = sample_info$factors,
    keep_num = 1,
    missing_value = c(0, NA)
  )

  sample_smd = sample_smd[keep_counts, ]
  sample_counts = sample_counts[keep_counts, ]

  if (id %in% "NSCLC") {
    sample_norm = sample_counts
  } else {
    sample_norm = median_normalize(sample_counts)
  }

  assays(sample_smd)$normalized = sample_norm

  sample_completeness = pairwise_completeness(sample_counts)
  ici_cor = ici_kendalltau(sample_norm, global_na = c(NA, 0))$cor

  sample_counts_na = sample_counts
  sample_counts_na[sample_counts_na == 0] = NA

  kt = kt_fast(
    sample_counts_na,
    use = "pairwise.complete.obs",
    return_matrix = TRUE
  )$tau
  # this one should match the Gierlinski paper values for median correlations
  pearson_base_nozero = cor(
    sample_counts_na,
    method = "pearson",
    use = "pairwise.complete.obs"
  )
  pearson_base = cor(
    sample_counts,
    method = "pearson",
    use = "pairwise.complete.obs"
  )
  pearson_log1p = cor(
    log1p(sample_counts),
    method = "pearson",
    use = "pairwise.complete.obs"
  )
  log_counts = log(sample_counts)
  log_counts[is.infinite(log_counts)] = NA
  pearson_log = cor(
    log_counts,
    method = "pearson",
    use = "pairwise.complete.obs"
  )

  cor_vals = list(
    icikt = ici_cor,
    icikt_complete = ici_cor * sample_completeness,
    pearson_base = pearson_base,
    pearson_base_nozero = pearson_base_nozero,
    pearson_log1p = pearson_log1p,
    pearson_log = pearson_log,
    kt_base = kt
  )
  out_data = list(
    id = id,
    data = sample_smd,
    ranks = rank_info,
    censorship = censorship_test,
    value_check = value_check,
    cor_vals = cor_vals
  )
  return(out_data)
}

run_kt = function(sample_counts) {
  sample_counts_na = sample_counts
  sample_counts_na[sample_counts_na == 0] = NA

  kt_fast(
    sample_counts_na,
    use = "pairwise.complete.obs",
    return_matrix = TRUE
  )$tau
}

run_cor_everyway = function(sample_counts, sample_completeness) {
  ici_cor = ici_kendalltau(sample_counts, global_na = c(NA, 0))$cor

  sample_norm_na = sample_norm
  sample_norm_na[sample_norm_na == 0] = NA

  kt = kt_fast(
    sample_norm_na,
    use = "pairwise.complete.obs",
    return_matrix = TRUE
  )$tau
  # this one should match the Gierlinski paper values for median correlations
  pearson_base_nozero = cor(
    sample_norm_na,
    method = "pearson",
    use = "pairwise.complete.obs"
  )
  pearson_base = cor(
    sample_norm,
    method = "pearson",
    use = "pairwise.complete.obs"
  )
  pearson_log1p = cor(
    log1p(sample_norm),
    method = "pearson",
    use = "pairwise.complete.obs"
  )
  log_counts = log(sample_norm)
  log_counts[is.infinite(log_counts)] = NA
  pearson_log = cor(
    log_counts,
    method = "pearson",
    use = "pairwise.complete.obs"
  )

  cor_vals = list(
    icikt = ici_cor,
    icikt_complete = ici_cor * sample_completeness,
    pearson_base = pearson_base,
    pearson_base_nozero = pearson_base_nozero,
    pearson_log1p = pearson_log1p,
    pearson_log = pearson_log,
    kt_base = kt
  )
  cor_vals
}

calculate_cor_medians = function(sample_cor, sample_ids, sample_classes) {
  out_med = purrr::imap(sample_cor, function(in_cor, cor_id) {
    intersect_ids = base::intersect(colnames(in_cor), sample_ids)
    keep_ids = sample_ids %in% intersect_ids
    use_ids = sample_ids[keep_ids]
    use_classes = sample_classes[keep_ids]

    in_cor = in_cor[use_ids, use_ids]
    in_med = visualizationQualityControl::median_correlations(
      in_cor,
      use_classes
    )
    in_med$which = cor_id
    in_med
  })
}


calculate_missingness = function(id, smd_file, type, value_check) {
  sample_smd = readRDS(smd_file)

  sample_counts = assays(sample_smd)$counts
  sample_info = colData(sample_smd) |> as.data.frame()

  keep_counts = keep_non_missing_percentage(
    sample_counts,
    sample_classes = sample_info$factors,
    keep_num = 1,
    missing_value = c(0, NA)
  )

  sample_smd = sample_smd[keep_counts, ]
  sample_counts = sample_counts[keep_counts, ]

  n_miss = sum(
    is.na(sample_counts) | (is.infinite(sample_counts)) | (sample_counts == 0)
  )
  n_value = nrow(sample_counts) * ncol(sample_counts)

  missing_table = tibble::tibble(
    id = id,
    type = type,
    missing = n_miss,
    values = n_value,
    value_check = value_check
  ) |>
    dplyr::mutate(perc_miss = missing / values * 100)
  missing_test = ICIKendallTau::test_left_censorship(
    sample_counts,
    sample_classes = sample_info$factors
  )
  missing_tidy = broom::tidy(missing_test$binomial_test)
  missing_table$p.value = missing_tidy$p.value
  missing_table
}

calculate_rank_correlation = function(id, smd_file, type, value_check) {
  sample_smd = readRDS(smd_file)

  sample_counts = assays(sample_smd)$counts
  sample_info = colData(sample_smd) |> as.data.frame()

  keep_counts = keep_non_missing_percentage(
    sample_counts,
    sample_classes = sample_info$factors,
    keep_num = 1,
    missing_value = c(0, NA)
  )

  sample_smd = sample_smd[keep_counts, ]
  sample_counts = sample_counts[keep_counts, ]

  ranked_data = ICIKendallTau::rank_order_data(
    sample_counts,
    sample_classes = sample_info$factors
  )

  rank_cor = purrr::imap(ranked_data, \(in_data, id) {
    # in_data = ranked_data[[1]]
    min_rank = in_data$n_na_rank |>
      dplyr::summarise(min_rank = min(median_rank), .by = n_na)
    full_cor = cor(
      in_data$n_na_rank$n_na,
      in_data$n_na_rank$median_rank,
      method = "kendall"
    )
    min_cor = cor(
      min_rank$n_na,
      min_rank$min_rank,
      method = "kendall"
    )
    tibble::tibble(median = full_cor, min = min_cor, factors = id)
  }) |>
    purrr::list_rbind()
  rank_cor$id = id
  rank_cor$type = type
  rank_cor$value_check = value_check

  rank_cor
}
