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
  row_rsd = data.frame(
    rsd = row_sd / row_mean,
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


calculate_outlier_effects = function(metabolomics_cor, type = "all") {
  # metabolomics_cor = tar_read(metabolomics_cor_AN000359)
  # metabolomics_cor = tar_read(metabolomics_cor_AN004368)
  # metabolomics_cor = tar_read(metabolomics_cor_AN001736)

  sample_counts = assays(metabolomics_cor$data)$normalized
  sample_info = colData(metabolomics_cor$data) |> as.data.frame()

  # instead of doing it this way, why not just do the whole set of
  # calculations twice in the workflow? Once with subsets,
  # and once without, using an argument to this function.
  if (type %in% "all") {
    raw_rsd = calculate_rsd_subset(sample_counts)
    raw_rsd$factors = "all"
    raw_rsd$correlation = "raw"
  } else if (type %in% "subsets") {
    raw_rsd = calculate_rsd(sample_counts, sample_info$factors)
    raw_rsd$correlation = "raw"
  }

  # these should be wrapped into functions, so they can be shorter calls, because
  # we want to do them twice.
  # 1 - figure out what features always have non-missing in 3 samples, regardless of what was removed,
  # 2 - subset to those features, and then do RSD calculations again.
  outlier_cor = purrr::imap(metabolomics_cor$cor, \(in_cor, cor_id) {
    # cor_id = "icikt"
    # in_cor = metabolomics_cor$cor[[cor_id]]
    # in_cor = metabolomics_cor$cor$icikt
    median_cor = median_correlations(in_cor, sample_info$factors)
    median_outlier = determine_outliers(
      median_correlations = median_cor,
      outlier_frac = NULL
    )
    median_outlier$correlation = cor_id
    median_outlier
  })
  outlier_rsd = purrr::imap(outlier_cor, \(in_outlier, cor_id) {
    # cor_id = "icikt"
    # in_outlier = outlier_cor[[cor_id]]

    keep_cor = in_outlier |> dplyr::filter(!outlier)
    keep_info = sample_info |> dplyr::filter(sample_id %in% keep_cor$sample_id)
    keep_counts = sample_counts[, keep_info$sample_id]

    if (type %in% "all") {
      keep_rsd = calculate_rsd_subset(keep_counts)
      keep_rsd$factors = "all"
    } else if (type %in% "subsets") {
      keep_rsd = calculate_rsd(keep_counts, sample_info$factors)
    }
    keep_rsd$correlation = cor_id
    keep_rsd
  }) |>
    purrr::list_rbind()
  outlier_rsd = dplyr::bind_rows(raw_rsd, outlier_rsd)

  metabolomics_cor$outliers = list(
    median_correlations = outlier_cor,
    rsd = outlier_rsd
  )
  return(metabolomics_cor)
}

calculate_rsd_differences = function(outliers) {
  # outliers = tar_read(outliers_all_NSCLC)
  outlier_rsd = outliers$outliers$rsd

  outlier_split = split(outlier_rsd, outlier_rsd$factors)

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
  })
  outliers$outliers$rsd_diffs = outlier_diff_raw
  return(outliers)
}

calculate_median_rsd_diffs = function(rsd_diffs) {
  # rsd_diffs = tar_read(rsd_all_AN003382)
  # rsd_diffs = tar_read(rsd_subsets_NSCLC)
  outlier_diff = rsd_diffs$outliers$rsd_diffs

  outlier_median = purrr::imap(outlier_diff, \(in_diff, id) {
    in_diff |>
      dplyr::summarise(
        diff_median = median(diff, na.rm = TRUE),
        perc_diff_median = median(perc_diff, na.rm = TRUE),
        .by = correlation
      ) |>
      dplyr::mutate(factors = id, id = rsd_diffs$id)
  }) |>
    purrr::list_rbind()
  rsd_diffs$outliers$rsd_diff_median = outlier_median
  rsd_diffs
}

get_just_medians = function(median_diffs) {
  # median_diffs = tar_read(rsd_all_median_AN003382)
  # median_diffs = tar_read(rsd_subsets_median_AN003382)
  median_diffs$outliers$rsd_diff_median
}
