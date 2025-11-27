create_lod_run_df = function(lod_combinations, lod_levels) {
  run_df = purrr::map(lod_combinations, \(in_comb) {
    # in_comb = lod_combinations[[3]]
    purrr::map(seq_len(ncol(in_comb)), \(in_col) {
      use_comb = in_comb[, in_col]
      tibble::tibble(
        multipliers = list(lod_levels$lod[use_comb]),
        id = paste0(lod_levels$level[use_comb], collapse = "_")
      )
    }) |>
      dplyr::bind_rows()
  }) |>
    dplyr::bind_rows()
  run_df
}

# an alternative method to limit of detection would be to use a varying dynamic
# range, where we vary from 3 - 2 orders of magnitude, and base it on the highest
# observed value. May be able to use a uniform random distribution.
#
# Also consider the differences not just as an absolute difference, but as a relative
# difference, where we take the difference / true correlation.

create_large_replicate_samples = function(lod_vars) {
  # tar_load(lod_vars)
  base_sample = rlnorm(
    lod_vars$n_feature,
    meanlog = lod_vars$meanlog,
    sdlog = lod_vars$sdlog
  )
  # should we increase meanlog = 3, sdlog = 1.5, then noise = 0.2
  rep_data = add_uniform_noise(lod_vars$n_sample, base_sample, lod_vars$sd)
  colnames(rep_data) = paste0(
    "S",
    stringr::str_pad(seq_len(ncol(rep_data)), width = 3, pad = "0")
  )

  log_data = log10(exp(rep_data))

  log_data
}


check_lod_na_perc = function(var_lod_samples, check_lod_levels) {
  # tar_load(var_lod_samples)
  out_na = purrr::map(colnames(var_lod_samples), \(in_sample) {
    sample_oom = min(var_lod_samples[, in_sample]) + check_lod_levels
    sample_na = purrr::map_dbl(sample_oom, \(in_oom) {
      (sum(var_lod_samples[, in_sample] < in_oom))
    })
    tibble::tibble(
      sample_id = in_sample,
      oom = check_lod_levels,
      sample_oom = sample_oom,
      sample_na = sample_na
    )
  }) |>
    purrr::list_rbind()

  out_na
}


create_variable_lod_samples = function(rep_data, oom_max, id) {
  # rep_data = tar_read(var_lod_samples)
  # oom_max = 1
  # id = "low"

  oom_vector = runif(ncol(rep_data), 0, oom_max)

  rep_mins = apply(rep_data, 2, min)
  lod_cutoffs = rep_mins + oom_vector

  lod_matrix = matrix(
    lod_cutoffs,
    nrow = nrow(rep_data),
    ncol = ncol(rep_data),
    byrow = TRUE
  )

  rep_lod = rep_data
  rep_lod[rep_data < lod_matrix] = NA

  n_na = apply(rep_lod, 2, \(x) {
    sum(is.na(x))
  })

  oom_info = tibble::tibble(
    sample_id = colnames(rep_data),
    min_value = rep_mins,
    oom = oom_vector,
    lod_cutoff = lod_cutoffs,
    n_na = n_na,
    oom_max = oom_max,
    oom_id = id
  )

  return(list(
    sample_data = rep_data,
    sample_lod = rep_lod,
    oom = oom_info,
    oom_id = id
  ))
}

variable_lod_cor_everyway = function(sample_counts) {
  # tmp_values = tar_read(var_lod_samples)
  # sample_counts = tmp_values$variable_cutoff
  # sample_counts = var_lod_samples$sample_data

  min_value = min(sample_counts, na.rm = TRUE)
  impute_value = min_value / 2

  sample_counts_min = sample_counts
  sample_counts_min[is.na(sample_counts)] = impute_value

  ici_cor = ici_kendalltau(
    sample_counts,
    global_na = NA,
    scale_max = FALSE
  )$cor
  kt_na = kt_fast(
    sample_counts,
    use = "pairwise.complete.obs",
    return_matrix = TRUE
  )$tau
  kt_min = kt_fast(
    sample_counts_min,
    use = "pairwise.complete.obs",
    return_matrix = TRUE
  )$tau
  pearson_na = cor(
    sample_counts,
    method = "pearson",
    use = "pairwise.complete.obs"
  )
  pearson_min = cor(
    sample_counts_min,
    method = "pearson",
    use = "pairwise.complete.obs"
  )

  cor_vals = list(
    icikt = ici_cor,
    kt_na = kt_na,
    kt_min = kt_min,
    pearson_na = pearson_na,
    pearson_min = pearson_min
  )
  cor_vals
}


calculate_variable_correlations = function(var_lod_samples) {
  # var_lod_samples = tar_read(vl_samples_low)
  # var_lod_samples = tar_read(vl_samples_high)
  impute_value = min(var_lod_samples$sample_lod, na.rm = TRUE) / 2
  ref_correlations = variable_lod_cor_everyway(var_lod_samples$sample_data)
  lod_correlations = variable_lod_cor_everyway(var_lod_samples$sample_lod)

  list(
    samples = var_lod_samples,
    oom = var_lod_samples$oom,
    reference_cor = ref_correlations,
    lod_cor = lod_correlations,
    impute_value = impute_value
  )
}

calculate_var_lod_correlation_diffs = function(var_lod_correlations) {
  # var_lod_correlations = tar_read(vl_cor_med)
  # var_lod_correlations = tar_read(vl_cor_0.5_1_1.25_1.5)
  reference_cor = var_lod_correlations$reference_cor
  lod_cor = var_lod_correlations$lod_cor

  diff_cor = purrr::map(names(reference_cor), \(cor_id) {
    diff_df = var_lod_each_cor(reference_cor[[cor_id]], lod_cor[[cor_id]])

    diff_df = diff_df |>
      dplyr::mutate(cor_method = cor_id)
    diff_df
  }) |>
    purrr::list_rbind()

  diff_cor_lod = add_lod_info(diff_cor, var_lod_correlations$oom)
  diff_cor_lod$impute_value = var_lod_correlations$impute_value
  diff_cor_lod$oom_id = var_lod_correlations$oom$oom_id[1]

  diff_cor_lod
}

add_lod_info = function(diff_cor, lod_df) {
  # lod_df = var_lod_correlations$oom

  diff_cor_lod = dplyr::left_join(
    diff_cor,
    lod_df |>
      dplyr::transmute(
        s1 = sample_id,
        s1_lod = lod_cutoff,
        s1_oom = oom,
        s1_min_value = min_value,
        s1_n_na = n_na
      ),
    by = "s1"
  )
  diff_cor_lod = dplyr::left_join(
    diff_cor_lod,
    lod_df |>
      dplyr::transmute(
        s2 = sample_id,
        s2_lod = lod_cutoff,
        s2_oom = oom,
        s2_min_value = min_value,
        s2_n_na = n_na
      ),
    by = "s2"
  )

  diff_cor_lod
}

var_lod_each_cor = function(ref_cor, in_cor) {
  # ref_cor = reference_cor[["icikt_na"]]
  # in_cor = var_lod_correlations[["single_cutoff"]][["icikt_na"]]
  ref_long = lod_cor_matrix_2_df(ref_cor) |>
    dplyr::mutate(comparison = glue::glue("{s1}.{s2}")) |>
    dplyr::filter(!(s1 == s2), !duplicated(comparison))
  in_long = lod_cor_matrix_2_df(in_cor) |>
    dplyr::mutate(comparison = glue::glue("{s1}.{s2}")) |>
    dplyr::filter(!(s1 == s2), !duplicated(comparison))

  compare_cor = dplyr::left_join(
    ref_long[, c("cor", "s1", "s2", "comparison")],
    in_long[, c("cor", "comparison")],
    suffix = c("_ref", "_lod"),
    by = "comparison"
  )
  compare_cor = compare_cor |>
    dplyr::mutate(
      ref_v_lod = cor_ref - cor_lod,
      ref_v_lod_relative = (cor_ref - cor_lod) / cor_ref
    )
  compare_cor
}

triple_check_scalemax_works = function(var_lod_samples) {
  use_sample = var_lod_samples$single_cutoff

  ici_noscale = ici_kendalltau(
    use_sample,
    global_na = NA,
    scale_max = FALSE,
    diag_good = FALSE
  )

  testthat::expect_equal(ici_noscale$cor, ici_noscale$raw)

  ici_scale = ici_kendalltau(
    use_sample,
    global_na = NA,
    scale_max = TRUE,
    diag_good = FALSE
  )
  sum_lt = sum(ici_scale$cor <= ici_scale$raw)
  testthat::expect_equal(sum_lt, nrow(ici_scale$cor) * ncol(ici_scale$cor))
  NULL
}

lod_cor_matrix_2_df = function(in_matrix) {
  comparisons = combn(colnames(in_matrix), 2)

  cor_vals = vector(mode = "double", length = ncol(comparisons))

  for (icol in seq_len(ncol(comparisons))) {
    cor_vals[icol] = in_matrix[comparisons[1, icol], comparisons[2, icol]]
  }

  cor_df = data.frame(
    s1 = comparisons[1, ],
    s2 = comparisons[2, ],
    cor = cor_vals
  )
  cor_df
}

lod_cor_sample_matrix_2_df = function(lod_cor) {
  # lod_cor = tar_read(vl_cor_all)
  sample_matrix = lod_cor$samples$sample_data
  sample_df = sample_matrix_2_df(sample_matrix)
  sample_df$lod_id = lod_cor$samples$sample_id
  sample_df = dplyr::left_join(sample_df, lod_cor$lod, by = "sample")
  sample_df$imputed_value = lod_cor$impute_value
  sample_df$which = "full"

  lod_matrix = lod_cor$samples$sample_lod
  lod_df = sample_matrix_2_df(lod_matrix)
  lod_df$lod_id = lod_cor$samples$sample_id
  lod_df = dplyr::left_join(lod_df, lod_cor$lod, by = "sample")
  lod_df$imputed_value = lod_cor$impute_value
  lod_df$which = "lod"

  dplyr::bind_rows(sample_df, lod_df)
}

sample_matrix_2_df = function(sample_matrix) {
  n_tot = nrow(sample_matrix) * ncol(sample_matrix)
  total_value = vector("numeric", n_tot)
  if (is.null(rownames(sample_matrix))) {
    row_names = paste0("f", seq_len(nrow(sample_matrix)))
  } else {
    row_names = rownames(sample_matrix)
  }

  row_id = vector("character", n_tot)
  col_id = vector("character", n_tot)
  i_entry = 1
  for (icol in colnames(sample_matrix)) {
    for (irow in seq_len(nrow(sample_matrix))) {
      row_id[i_entry] = row_names[irow]
      col_id[i_entry] = icol
      total_value[i_entry] = sample_matrix[irow, icol]
      i_entry = i_entry + 1
    }
  }
  tibble::tibble(feature = row_id, sample = col_id, value = total_value)
}

calculate_cor_diff_summaries = function(vl_cor_diff) {
  # vl_cor_diff = tar_read(vl_cor_diff_0.5_1.5)
  vl_summary = vl_cor_diff |>
    dplyr::group_by(cor_method, compare_levels) |>
    dplyr::summarise(
      ref_v_lod_median = median(abs(ref_v_lod)),
      s1_percna_median = median(s1_perc_na),
      s2_percna_median = median(s2_perc_na),
      n_lod = n_lod[1],
      lod_id = lod_id[1]
    )
  vl_summary
}


create_na_perc_graph = function(vl_na_perc, lod_ranges_tar) {
  # tar_load(c(vl_na_perc,
  #             lod_ranges_tar))
  vl_na_perc |>
    ggplot(aes(x = oom, y = sample_na, group = oom)) +
    geom_sina() +
    geom_vline(xintercept = c(lod_ranges_tar$max), color = "red") +
    labs(x = "Change in Dynamic Range", y = "# Missing Values")
}

create_lod_diff_graph = function(vl_cor_diff_all) {
  # tar_load(vl_cor_diff_all)
  icikt_min = vl_cor_diff_all |>
    dplyr::filter(cor_method %in% c("icikt", "pearson_min")) |>
    dplyr::mutate(oom_id = factor(oom_id, levels = c("low", "med", "high")))

  icikt_min_wider = icikt_min |>
    dplyr::select(
      oom_id,
      comparison,
      ref_v_lod,
      s1_n_na,
      s2_n_na,
      cor_method
    ) |>
    tidyr::pivot_wider(names_from = cor_method, values_from = ref_v_lod) |>
    dplyr::mutate(
      max_na = dplyr::case_when(
        s1_n_na >= s2_n_na ~ s1_n_na,
        s1_n_na <= s2_n_na ~ s2_n_na
      )
    )

  # split into each level, and do it via patchwork so each one can have it's own
  # coordinates.
  split_level = split(icikt_min_wider, icikt_min_wider$oom_id)
  full_na_counts = seq(1, max(icikt_min_wider$max_na), 1)
  na_count_range = range(full_na_counts)
  na_count_scale = scales::rescale(full_na_counts)
  against_graphs = purrr::map(split_level, \(in_level) {
    if (in_level$oom_id[1] %in% "high") {
      show_legend = TRUE
      out_graph = in_level |>
        ggplot(aes(x = icikt, y = pearson_min, color = max_na)) +
        geom_abline(slope = 1, color = "red") +
        geom_abline(slope = -1, color = "red") +
        geom_point(alpha = 0.25, show.legend = show_legend) +
        scale_color_viridis_c(
          values = na_count_scale,
          limits = na_count_range
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        labs(
          x = "ICI-Kt (ref - trimmed)",
          y = "Pearson Imputed (ref - trimmed)",
          color = "# NA",
          subtitle = stringr::str_to_title(in_level$oom_id[1])
        ) +
        theme(legend.position = c(0.85, 0.22))
    } else {
      show_legend = FALSE
      out_graph = in_level |>
        ggplot(aes(x = icikt, y = pearson_min, color = max_na)) +
        geom_abline(slope = 1, color = "red") +
        geom_abline(slope = -1, color = "red") +
        geom_point(alpha = 0.25, show.legend = show_legend) +
        scale_color_viridis_c(
          values = na_count_scale,
          limits = na_count_range
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        labs(
          x = "ICI-Kt (ref - trimmed)",
          y = "Pearson Imputed (ref - trimmed)",
          subtitle = stringr::str_to_title(in_level$oom_id[1])
        )
    }
    out_graph
  })

  hist_graphs = purrr::map(split_level, \(in_level) {
    in_level |>
      ggplot(aes(x = abs(icikt) - abs(pearson_min))) +
      geom_histogram(bins = 100) +
      geom_vline(xintercept = 0, color = "red") +
      labs(
        x = "abs(ICI-Kt) - abs(Pearson Min)",
        subtitle = stringr::str_to_title(in_level$oom_id[1])
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, .1)))
  })

  diff_hist = c(against_graphs, hist_graphs)
  diff_hist
}

create_icikt_ktimpute_graph = function(vl_cor_diff_all) {
  # tar_load(vl_cor_diff_all)
  icikt_min = vl_cor_diff_all |>
    dplyr::filter(cor_method %in% c("icikt", "kt_min"), oom_id %in% "high")

  icikt_min_wider = icikt_min |>
    dplyr::select(oom_id, comparison, cor_lod, cor_method) |>
    tidyr::pivot_wider(names_from = cor_method, values_from = cor_lod)

  icikt_min_wider |>
    ggplot(aes(x = icikt, y = kt_min)) +
    geom_abline(slope = 1, color = "red") +
    geom_point() +
    labs(x = "ICI-Kt", y = "Kt Imputed")
}

create_lod_hist_graph = function(vl_cor_diff_all) {
  # tar_load(vl_cor_diff_all)
  NULL
}
