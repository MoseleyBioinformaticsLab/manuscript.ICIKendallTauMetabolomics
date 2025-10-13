filter_generate_outliers = function(
  counts,
  info,
  keep_num,
  sample_col,
  class_col
) {
  # counts_info = tar_read(yeast_counts_info)
  # counts = counts_info$counts
  # info = counts_info$info
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "treatment"
  #
  #
  # counts_info = tar_read(nsclc_counts_info)
  # counts = counts_info$counts
  # info = counts_info$info
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "treatment"
  if (length(class_col) == 2) {
    filter_col = class_col[1]
    median_col = class_col[2]
  } else {
    filter_col = class_col
    median_col = class_col
  }
  counts_filter = t(keep_non_zero_percentage(
    t(counts),
    sample_classes = info[[filter_col]],
    keep_num = keep_num
  ))
  counts_completeness = pairwise_completeness(t(counts_filter))
  counts_cor = run_cor_everyway(counts_filter, counts_completeness)
  counts_medians = calculate_cor_medians(
    counts_cor,
    info[[sample_col]],
    info[[median_col]]
  )

  counts_outliers = purrr::map_df(counts_medians, determine_outliers)
  counts_outliers$keep_num = keep_num
  list(outliers = counts_outliers, features = rownames(counts))
}


get_single_outlier = function(outlier_list) {
  # outlier_list = tar_read(yeast_outliers)
  n_frac = purrr::map_dbl(outlier_list, \(x) {
    x$outliers$keep_num[1]
  })
  outlier_list[[which(n_frac == 1)]]$outliers
}

rbind_outliers = function(outlier_list) {
  all_outliers = purrr::map(outlier_list, \(x) {
    x$outliers
  }) |>
    dplyr::bind_rows()
  all_outliers
}

calculate_rsd_subset = function(in_matrix) {
  row_mean = base::rowMeans(in_matrix, na.rm = TRUE)
  row_sd = apply(in_matrix, 1, sd, na.rm = TRUE)
  row_n = apply(in_matrix, 1, \(x) {
    sum(!is.na(x))
  })
  row_rsd = data.frame(
    rsd = row_sd / row_mean,
    n = row_n,
    feature_id = rownames(in_matrix)
  )
  row_rsd
}

calculate_rsd = function(in_counts, in_classes) {
  # in_counts = sample_counts
  # in_classes = sample_info$factors

  split_samples = split(colnames(in_counts), in_classes)

  rsd_values = purrr::imap(split_samples, \(in_samples, split_id) {
    use_subset = in_counts[, in_samples, drop = FALSE]
    if (ncol(use_subset) < 3) {
      return(NULL)
    }
    rsd_subset = calculate_rsd_subset(use_subset)

    rsd_subset$factors = split_id
    rsd_subset
  }) |>
    purrr::list_rbind()

  rsd_values
}

calculate_raw_rsd = function(sample_counts, sample_info, type = "all") {
  if (type %in% "all") {
    raw_rsd = calculate_rsd_subset(sample_counts)
    raw_rsd$factors = "all"
    raw_rsd$correlation = "raw"
  } else if (type %in% "subsets") {
    raw_rsd = calculate_rsd(sample_counts, sample_info$factors)
    raw_rsd$correlation = "raw"
  }
  raw_rsd
}

calculate_outlier_rsd = function(
  outlier_cor,
  sample_counts,
  sample_info,
  type = "all"
) {
  # outlier_cor = metabolomics_cor$cor_vals
  found_outliers = purrr::imap(outlier_cor, \(in_cor, cor_id) {
    # cor_id = "icikt"
    # in_cor = outlier_cor[[cor_id]]
    # in_cor = metabolomics_cor$cor$icikt
    median_cor = median_correlations(in_cor, sample_info$factors)
    median_outlier = determine_outliers(
      median_correlations = median_cor,
      outlier_frac = NULL
    )
    median_outlier$correlation = cor_id
    median_outlier
  })

  outlier_rsd = purrr::imap(found_outliers, \(in_outlier, cor_id) {
    # cor_id = "icikt"
    # in_outlier = found_outliers[[cor_id]]

    keep_cor = in_outlier |> dplyr::filter(!outlier)
    keep_info = sample_info |> dplyr::filter(sample_id %in% keep_cor$sample_id)
    keep_counts = sample_counts[, keep_info$sample_id]

    if (type %in% "all") {
      keep_rsd = calculate_rsd_subset(keep_counts)
      keep_rsd$factors = "all"
    } else if (type %in% "subsets") {
      keep_rsd = calculate_rsd(keep_counts, keep_info$factors)
    }
    keep_rsd$correlation = cor_id
    keep_rsd
  }) |>
    purrr::list_rbind()
  outlier_rsd
}


calculate_outlier_effects = function(
  metabolomics_cor,
  metabolomics_keep,
  type = "all"
) {
  # metabolomics_cor = tar_read(metabolomics_cor_AN000359)
  # metabolomics_cor = tar_read(metabolomics_cor_AN004368)
  # metabolomics_cor = tar_read(metabolomics_cor_AN001736)

  # metabolomics_cor = tar_read(metabolomics_cor_AN001156)
  # metabolomics_keep = tar_read(metabolomics_keep_AN001156)
  # type = "all"
  sample_counts = assays(metabolomics_keep)$normalized
  sample_info = colData(metabolomics_keep) |> as.data.frame()

  raw_rsd = calculate_raw_rsd(sample_counts, sample_info, type)

  # these should be wrapped into functions, so they can be shorter calls, because
  # we want to do them twice.
  # 1 - figure out what features always have non-missing in 3 samples, regardless of what was removed,
  # 2 - subset to those features, and then do RSD calculations again.

  outlier_rsd = calculate_outlier_rsd(
    metabolomics_cor$cor_vals,
    sample_counts,
    sample_info,
    type = type
  )

  n_factors = length(unique(outlier_rsd$factors))
  n_cor = length(unique(outlier_rsd$correlation))
  outlier_n = outlier_rsd |>
    dplyr::summarise(
      has_n = sum(n >= 3) == n_factors,
      .by = c(feature_id, correlation)
    ) |>
    dplyr::summarise(
      has_cor = sum(has_n) == n_cor,
      .by = feature_id
    )

  tmp_n = outlier_n |>
    dplyr::filter(has_cor)

  if (nrow(tmp_n) == 0) {
    return(NULL)
  }

  raw_rsd = raw_rsd |>
    dplyr::filter(feature_id %in% tmp_n$feature_id)
  outlier_rsd = outlier_rsd |>
    dplyr::filter(feature_id %in% tmp_n$feature_id)

  all_rsd = dplyr::bind_rows(raw_rsd, outlier_rsd)
  return(all_rsd)
}

calculate_rsd_differences = function(outliers) {
  # outliers = tar_read(outliers_all_NSCLC)
  # outliers = tar_read(outliers_all_AN002987)
  # outliers = tar_read(outliers_subsets_AN004314)
  if (is.null(outliers)) {
    return(NULL)
  }
  outlier_split = split(outliers, outliers$factors)

  outlier_diff_raw = purrr::map(outlier_split, \(in_split) {
    # in_split = outlier_split[[1]]
    in_raw = in_split |>
      dplyr::filter(correlation %in% "raw") |>
      dplyr::mutate(other = rsd)

    diff_rsd = dplyr::left_join(
      in_split,
      in_raw[, c("feature_id", "other")],
      by = "feature_id"
    )

    diff_rsd = diff_rsd |>
      dplyr::mutate(diff = other - rsd, perc_diff = diff / other)
    diff_rsd
  }) |>
    purrr::list_rbind()
  return(outlier_diff_raw)
}

calculate_median_rsd_diffs = function(rsd_diffs) {
  # rsd_diffs = tar_read(rsd_all_AN003382)
  # rsd_diffs = tar_read(rsd_subsets_NSCLC)
  if (is.null(rsd_diffs)) {
    return(NULL)
  }

  diff_median = rsd_diffs |>
    dplyr::summarize(
      diff_median = median(diff, na.rm = TRUE),
      diff_mad = mad(diff, na.rm = TRUE),
      perc_diff_median = median(perc_diff, na.rm = TRUE),
      perc_diff_mad = mad(perc_diff, na.rm = TRUE),
      .by = c(correlation, factors)
    )

  diff_median
}

get_just_medians = function(median_diffs) {
  # median_diffs = tar_read(rsd_all_median_AN003382)
  # median_diffs = tar_read(rsd_subsets_median_AN003382)
  median_diffs$outliers$rsd_diff_median
}

filter_outlier_dopca = function(metabolomics_cor, metabolomics_keep) {
  # metabolomics_cor = tar_read(metabolomics_cor_AN001156)
  # metabolomics_keep = tar_read(metabolomics_keep_AN001156)

  sample_counts = assays(metabolomics_keep)$normalized
  sample_info = colData(metabolomics_keep) |> as.data.frame()

  min_counts = min(sample_counts, na.rm = TRUE) / 2
  sample_imputed = sample_counts
  sample_imputed[is.na(sample_counts)] = min_counts

  sample_log = log2(sample_imputed)

  keep_samples = purrr::map(
    metabolomics_cor$cor_vals,
    find_remove_outliers,
    sample_info
  )

  # don't forget to check without removing any outliers
  keep_samples$original = sample_info$sample_id

  n_org = nrow(sample_info)
  n_remove = purrr::imap(keep_samples, \(x, id) {
    tibble::tibble(correlation = id, n_removed = n_org - length(x))
  }) |>
    purrr::list_rbind()

  consider_pcs = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")
  pca_anova = purrr::imap(keep_samples, \(in_samples, cor_id) {
    dopca_test(in_samples, sample_log, sample_info, cor_id, consider_pcs)
  }) |>
    purrr::list_rbind()

  pc_likely = pca_anova |>
    dplyr::filter(correlation %in% "original", PC %in% consider_pcs) |>
    dplyr::arrange(p.value) |>
    dplyr::slice_head(n = 1) |>
    dplyr::pull(PC)

  pc_anova_imp = pca_anova |>
    dplyr::filter(PC %in% pc_likely)

  pc_anova_imp = cbind(pc_anova_imp, metadata(metabolomics_keep)$other)

  return(pc_anova_imp)
}

find_remove_outliers = function(in_cor, sample_info) {
  # in_cor = metabolomics_cor$cor_vals[[1]]
  in_cor = in_cor[sample_info$sample_id, sample_info$sample_id]

  med_cor = visualizationQualityControl::median_correlations(
    in_cor,
    sample_info$factors
  )
  outlier_status = visualizationQualityControl::determine_outliers(
    med_cor,
    outlier_fraction = NULL
  )

  keep_samples = outlier_status |>
    dplyr::filter(!outlier) |>
    dplyr::pull(sample_id)
  keep_samples
}


dopca_test = function(
  keep_samples,
  sample_log,
  sample_info,
  correlation,
  consider_pcs
) {
  sample_info = sample_info |>
    dplyr::filter(sample_id %in% keep_samples)
  sample_log = sample_log[, sample_info$sample_id]

  sample_pca = prcomp(t(sample_log), center = TRUE, scale. = FALSE)
  sample_scores = sample_pca$x[, consider_pcs]

  anova_res = custom_test_pca_scores(
    sample_scores,
    sample_info[, c("factors"), drop = FALSE]
  )
  anova_res$correlation = correlation
  anova_res
}


custom_test_pca_scores = function(pca_scores, sample_info) {
  pc_test = colnames(pca_scores)
  pc_2_variable = purrr::imap_dfr(sample_info, function(var_col, var_name) {
    pc_test = purrr::map_df(pc_test, function(in_pc) {
      #message(paste0(var_name, " : ", in_pc))
      tmp_col = var_col
      tmp_col = tmp_col[
        !is.na(var_col) & !is.infinite(var_col) & !is.nan(var_col)
      ]
      n_var = length(unique(tmp_col))
      is_1_all = FALSE
      if (n_var == 1) {
        is_1_all = TRUE
      }
      if ((n_var == nrow(sample_info)) && (is.character(var_col))) {
        is_1_all = TRUE
      }

      if (!is_1_all) {
        tmp_frame = data.frame(y = pca_scores[, in_pc], x = var_col)
        na_x = is.na(tmp_frame$x) |
          is.infinite(tmp_frame$x) |
          is.nan(tmp_frame$x)
        tmp_frame = tmp_frame[!na_x, ]
        aov_res = stats::aov(y ~ x, data = tmp_frame)
        tidy_res = broom::tidy(aov_res)[1, ]
        eta_val = effectsize::eta_squared(aov_res, partial = FALSE)
        tidy_res$eta2 = eta_val$Eta2
        tidy_res$PC = in_pc
        tidy_res$variable = var_name
        return(tidy_res)
      } else {
        return(NULL)
      }
    })
  })
  pc_2_variable
}


pca_compare_original = function(pca_outliers_all) {
  # tar_load(pca_outliers_all)
  keep_eta2 = pca_outliers_all |>
    dplyr::summarize(n_eta2 = length(unique(eta2)), .by = id) |>
    dplyr::filter(n_eta2 > 1)
  out_org = pca_outliers_all |>
    dplyr::filter(id %in% keep_eta2$id) |>
    dplyr::filter(correlation %in% "original") |>
    dplyr::select(id, correlation, eta2) |>
    dplyr::mutate(et2_org = eta2) |>
    dplyr::select(-eta2, -correlation)

  pca_eta2 = dplyr::inner_join(
    pca_outliers_all |>
      dplyr::select(id, eta2, correlation),
    out_org,
    by = "id"
  )

  pca_eta2 = pca_eta2 |>
    dplyr::mutate(eta2_diff = eta2 - et2_org)

  split_eta2 = split(pca_eta2, pca_eta2$id)

  pca_eta_rank = purrr::map(split_eta2, \(in_eta) {
    #in_eta = split_eta2[[1]]
    rank_vals = tibble::tibble(
      eta2 = sort(unique(in_eta$eta2), decreasing = TRUE)
    )
    rank_vals$rank = seq_len(nrow(rank_vals))
    out_eta = dplyr::left_join(in_eta, rank_vals, by = "eta2")
    out_eta
  }) |>
    purrr::list_rbind()

  rank_figure = pca_eta_rank |>
    ggplot(aes(x = rank, fill = correlation)) +
    geom_histogram(position = position_dodge(width = 0.5)) +
    theme(legend.position = "inside", legend.position.inside = c(0.7, 0.7)) +
    labs(x = TeX("$rank(\\eta^2)$"))

  full_histogram = pca_eta_rank |>
    dplyr::filter(!(correlation %in% "original")) |>
    ggplot(aes(y = eta2_diff, fill = correlation)) +
    geom_histogram(bins = 100) +
    facet_wrap(~correlation, nrow = 1) +
    theme(legend.position = "none") +
    labs(y = TeX("$\\eta^2 - \\eta{^2}_{org}$"))

  zoom_histogram = full_histogram +
    coord_cartesian(xlim = c(0, 20))

  list(rank = rank_figure, full = full_histogram, zoom = zoom_histogram)
}
