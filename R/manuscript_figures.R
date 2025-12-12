create_missing_percentage_sina_plot = function(missingness_summary) {
  # tar_load(missingness_summary)
  out_plot = missingness_summary |>
    ggplot(aes(x = perc_miss, y = type)) +
    geom_sina(size = 2) +
    labs(x = "Percent Missingness", y = "Metabolomics Method")
  return(out_plot)
}

create_rank_missingness_relationship_plot = function(in_ranks_metadata) {
  # in_cor = tar_read("metabolomics_cor_NSCLC")
  # in_ranks_metadata = tar_read("missingness_ranks_AN001074")
  rank_counts = in_ranks_metadata$ranked_data

  rank_counts_df = purrr::map(rank_counts, \(x) {
    x$n_na_rank
  }) |>
    purrr::list_rbind()
  rank_counts_min = rank_counts_df |>
    dplyr::summarise(min_rank = min(median_rank), .by = c(split, n_na))
  rank_na_image = rank_counts_df |>
    ggplot(aes(x = n_na, y = median_rank)) +
    geom_point(size = 2) +
    geom_point(
      data = rank_counts_min,
      aes(x = n_na, y = min_rank),
      color = "red"
    ) +
    facet_wrap(~split) +
    labs(x = "Number of Samples Missing In", y = "Median Rank") +
    theme(strip.text = element_text(size = 10))

  rank_na_image
}


create_rank_missingness_correlation_plot = function(rank_correlations) {
  # tar_load(rank_correlations)

  out_plot = rank_correlations |>
    tidyr::pivot_longer(
      c(median, min),
      values_to = "correlation",
      names_to = "summary_type"
    ) |>
    dplyr::mutate(
      rank_type = case_match(
        summary_type,
        "median" ~ "Median",
        "min" ~ "Min(median)"
      ),
      rank_type = factor(rank_type, levels = c("Min(median)", "Median")),
      type = factor(type, levels = c("NMR", "MS"))
    ) |>
    ggplot(aes(x = correlation, y = rank_type)) +
    geom_sina(size = 2, alpha = 0.5) +
    facet_wrap(~type, ncol = 1) +
    labs(x = "Kendall-tau(Rank vs N-Missing)", y = "Rank Summary")

  out_plot
}

create_missingness_pvalue_plot = function(missingness_summary) {
  # tar_load(missingness_summary)
  min_nonzero = missingness_summary |>
    dplyr::filter(padjust > 0) |>
    dplyr::pull(padjust) |>
    min()
  missingness_summary = missingness_summary |>
    dplyr::mutate(
      padjust2 = dplyr::case_when(
        padjust == 0 ~ min_nonzero,
        TRUE ~ padjust
      )
    )
  p_cut = -1 * log10(0.05)
  out_plot = missingness_summary |>
    dplyr::mutate(log_p = -1 * log10(padjust2)) |>
    ggplot(aes(x = log_p)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = p_cut, color = "red") +
    labs(x = "-1xLog10(P-adjusted)", y = "# Datasets") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0))
  out_plot
}

basic_ragg = function(
  filename,
  width = 8,
  height = 6,
  units = "in",
  res = 600,
  ...
) {
  ragg::agg_png(
    filename = filename,
    width = width,
    height = height,
    units = units,
    res = res,
    ...
  )
}

create_variable_dynamic_range_image = function(vl_diff_graph) {
  basic_ragg(
    "docs/poster/images/variable_dynamic_range.png",
    height = 8,
    width = 16
  )
  print(wrap_plots(vl_diff_graph, nrow = 2, ncol = 3))
  dev.off()
  return(invisible(NULL))
}

create_censored_value_images = function(censored_value_plots) {
  censored_value_plots = purrr::map(censored_value_plots, \(x) {
    x +
      theme(
        strip.background = element_rect(fill = "white", color = "white"),
        panel.border = element_blank()
      )
  })

  basic_ragg(
    "docs/poster/images/censored_value_plot.png",
    width = 16,
    height = 8
  )

  print(
    wrap_plots(censored_value_plots, ncol = 1) +
      plot_annotation(tag_levels = "A")
  )
  dev.off()
  return(invisible(NULL))
}

create_pca_outlier_image = function(pca_outlier_comparisons) {
  # tar_load(pca_outlier_comparisons)
  basic_ragg(
    "docs/poster/images/pca_outlier_plot.png",
    width = 16,
    height = 14
  )
  print(
    wrap_plots(pca_outlier_comparisons, ncol = 1) +
      plot_annotation(tag_levels = "A")
  )
  dev.off()
  return(invisible(NULL))
}


create_figure1_image = function(
  missingness_summary,
  example_rank,
  rank_correlations
) {
  # tar_load(missingness_summary)
  # example_rank = tar_read(missingness_tests_AN001074)
  # tar_load(rank_correlations)

  missingness_plot = missingness_summary |>
    ggplot(aes(x = perc_miss, y = type)) +
    geom_sina(size = 2) +
    labs(x = "% Missing Value", y = "Method")

  rank_cor_plot = rank_correlations |>
    tidyr::pivot_longer(
      c(median, min),
      values_to = "correlation",
      names_to = "summary_type"
    ) |>
    dplyr::mutate(
      rank_type = case_match(
        summary_type,
        "median" ~ "Median",
        "min" ~ "Min(median)"
      )
    ) |>
    ggplot(aes(x = correlation, y = rank_type)) +
    geom_sina(size = 2) +
    facet_wrap(~type, ncol = 1) +
    labs(x = "Kendall-tau(Rank vs N-Missing)", y = "Rank Summary")

  use_ranks = example_rank$rank_data |>
    purrr::map(\(x) {
      x$n_na_rank
    }) |>
    purrr::list_rbind() |>
    dplyr::mutate(factors = split)

  use_ranks_min = use_ranks |>
    dplyr::summarise(median_rank = min(median_rank), .by = c(factors, n_na))

  median_image = use_ranks |>
    ggplot(aes(x = n_na, y = median_rank)) +
    geom_point(size = 2, alpha = 0.5) +
    geom_point(data = use_ranks_min, size = 2, color = "red") +
    facet_wrap(~factors) +
    labs(x = "Number of Samples Missing In", y = "Median Rank") +
    theme(strip.text = element_text(size = 10))

  out_image = (missingness_plot | median_image) /
    rank_cor_plot +
    plot_annotation(tag_levels = "A")

  basic_ragg(
    "docs/poster/images/missingness_stuff.png",
    width = 16,
    height = 10
  )
  print(out_image)
  dev.off()
  invisible(NULL)
}

create_tabular_comparison_outputs = function(limma_comparisons) {
  # tar_load(limma_comparisons)
  frac_total = limma_comparisons$fraction_total |>
    dplyr::mutate(overall = sig_frac)

  order_values = frac_total |>
    dplyr::filter(overall %in% "high") |>
    dplyr::arrange(dplyr::desc(mean)) |>
    dplyr::select(method)

  high_values = dplyr::left_join(
    order_values,
    frac_total |> dplyr::filter(overall %in% "high"),
    by = "method"
  )
  high_gt = high_values |>
    dplyr::select(method, mean, sd) |>
    dplyr::transmute(Method = method, Mean = mean, SD = sd) |>
    gt::gt() |>
    gt::fmt_number(decimals = 4) |>
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        rows = c(1, 3)
      )
    )

  low_values = dplyr::left_join(
    order_values,
    frac_total |> dplyr::filter(overall %in% "low"),
    by = "method"
  )
  low_gt = low_values |>
    dplyr::select(method, mean, sd) |>
    dplyr::transmute(Method = method, Mean = mean, SD = sd) |>
    gt::gt() |>
    gt::fmt_number(decimals = 4) |>
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        rows = c(1, 3)
      )
    )

  # high_gtable = wrap_elements(gt::as_gtable(high_gt))
  # low_gtable = wrap_elements(gt::as_gtable(low_gt))

  # out_image = (high_gtable | low_gtable) + plot_annotation(tag_levels = "A")

  gt::gtsave(
    high_gt,
    filename = "docs/poster/images/high_stats.png",
    zoom = 2
  )
  gt::gtsave(low_gt, filename = "docs/poster/images/low_stats.png", zoom = 2)
  return(invisible(NULL))
}

create_tabular_test_outputs = function(test_comparisons) {
  # test_comparisons = tar_read(limma_test_comparisons_good)

  table_comparisons = test_comparisons |>
    dplyr::filter(p.value <= 0.05) |>
    dplyr::arrange(p.value) |>
    dplyr::transmute(
      Comparison = gsub("_v_", "\nvs\n", comparison),
      `P-Value` = p.value,
      Difference = estimate
    )
  table_gt = table_comparisons |>
    gt::gt() |>
    gt::fmt_scientific(decimals = 1) |>
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        rows = c(1, 2)
      )
    )

  gt::gtsave(table_gt, file = "docs/poster/images/t_test_stats.png")
  return(invisible(NULL))
}

create_graphical_comparison_outputs = function(limma_comparisons) {
  # tar_load(limma_comparisons)
  frac_total = limma_comparisons$fraction_total |>
    dplyr::mutate(overall = stringr::str_to_sentence(sig_frac))

  order_values = frac_total |>
    dplyr::filter(overall %in% "High") |>
    dplyr::arrange((mean))

  frac_total = frac_total |>
    dplyr::mutate(low = mean - (sd), high = mean + (sd))
  frac_total$method = factor(frac_total$method, levels = order_values$method)
  frac_total$overall = factor(frac_total$overall, levels = c("Low", "High"))
  mean_plot = frac_total |>
    ggplot(aes(x = mean, y = method)) +
    geom_point(size = 3) +
    facet_wrap(~overall, nrow = 1, scales = "free_x") +
    labs(x = "Mean(Fraction Significant)", y = "Method")

  sd_plot = frac_total |>
    ggplot(aes(x = sd, y = method)) +
    geom_point(size = 3) +
    facet_wrap(~overall, nrow = 1, scales = "free_x") +
    labs(x = "SD(Fraction Significant)", y = "Method")

  basic_ragg(filename = "docs/poster/images/mean_comparison.png")
  print(mean_plot)
  dev.off()

  basic_ragg(filename = "docs/poster/images/sd_comparison.png")
  print(sd_plot)
  dev.off()

  return(invisible(NULL))
}


create_graphical_test_outputs = function(test_comparisons) {
  # test_comparisons = tar_read(limma_test_comparisons_good)

  table_comparisons = test_comparisons |>
    dplyr::filter(p.value <= 0.05) |>
    dplyr::arrange(dplyr::desc(p.value)) |>
    dplyr::transmute(
      Comparison = gsub("_v_", "\nvs\n", comparison),
      `P-Value` = p.value,
      Difference = estimate,
      Low = conf.low,
      High = conf.high,
      logP = -1 * log10(`P-Value`)
    )

  table_comparisons$Comparison = factor(
    table_comparisons$Comparison,
    levels = table_comparisons$Comparison
  )

  table_plot = table_comparisons |>
    ggplot(aes(x = Difference, y = Comparison)) +
    geom_point(aes(size = logP)) +
    geom_linerange(aes(xmin = Low, xmax = High)) +
    theme(legend.position = "none")

  basic_ragg("docs/poster/images/t_test_graphic.png")
  print(table_plot)
  dev.off()
  return(invisible(NULL))
}


create_filtering_flowchart = function(check_data) {
  # out_graph = DiagrammeR::mermaid(
  #   "
  # flowchart TB
  # A[Repaired Datasets\nN = 6105] --> B(Have Metabolites\nN = 4469)
  # B --> C(N Metabolites >= 100 \n N = 2307)
  # C --> D(N Samples >= 5\n N SSF >= 2\n N = 1552)
  # D --> E(Max Intensity >= 20\nN = 1480)
  # E --> F[Correlation of N-Missing vs Rank != NA\nN = 711]
  # "
  # )

  check_fc = flowchart::as_fc(check_data, label = "Repaired Datasets") |>
    flowchart::fc_filter(
      !(FEATURE_CHECK %in% "NA FEATURES"),
      label = "Have Metabolite Data",
      show_exc = TRUE,
      perc_total = TRUE
    ) |>
    flowchart::fc_filter(
      FEATURE_CHECK %in% "GOOD",
      label = "N Metabolites >= 100",
      show_exc = TRUE,
      perc_total = TRUE
    ) |>
    flowchart::fc_filter(
      SSF_CHECK %in% "GOOD SSF",
      label = "N Samples >= 5\nN SSF >= 2",
      show_exc = TRUE,
      perc_total = TRUE
    ) |>
    flowchart::fc_filter(
      RANGE_CHECK %in% "GOOD",
      label = "Max Intensity >= 20",
      show_exc = TRUE,
      perc_total = TRUE
    ) |>
    flowchart::fc_filter(
      RANK_CHECK %in% c("SIGN DIFFERENCE", "GOOD"),
      label = "Correlation of\nN-Missing vs Rank != NA",
      show_exc = TRUE,
      perc_total = TRUE
    )
  check_fc |>
    flowchart::fc_draw() |>
    flowchart::fc_export(
      filename = here::here("docs/filtering_flowchart.png"),
      format = "png",
      width = 4,
      height = 10,
      units = "in",
      res = 300
    )
  return_file_hash(here::here("docs/filtering_flowchart.png"))
}

create_missingness_percentage_pvalue_plot = function(missingness_summary) {
  # tar_load(missingness_summary)
  min_nonzero = missingness_summary |>
    dplyr::filter(padjust > 0) |>
    dplyr::pull(padjust) |>
    min()
  missingness_summary = missingness_summary |>
    dplyr::mutate(
      padjust2 = dplyr::case_when(
        padjust == 0 ~ min_nonzero,
        TRUE ~ padjust
      )
    ) |>
    dplyr::mutate(log_p = -1 * log10(padjust2))

  out_plot = missingness_summary |>
    ggplot(aes(x = perc_miss, y = log_p)) +
    geom_point(size = 2, alpha = 0.5) +
    labs(x = "Percent Missingness", y = "-1xLog10(P-adjusted)")
  out_plot
}

return_file_hash = function(file_loc) {
  digest::digest(file_loc, algo = "sha256", file = TRUE)
}
count_filtered_datasets = function(check_data) {
  n_total = nrow(check_data)
  has_features = check_data |>
    dplyr::filter(!(FEATURE_CHECK %in% "NA METABOLITES"))
  good_features = check_data |>
    dplyr::filter(FEATURE_CHECK %in% "GOOD")
  good_ssf = good_features |>
    dplyr::filter(SSF_CHECK %in% "GOOD SSF")
  good_range = good_ssf |>
    dplyr::filter(RANGE_CHECK %in% "GOOD")
  good_rank = good_range |>
    dplyr::filter(RANK_CHECK %in% c("SIGN DIFFERENCE", "GOOD"))

  n_filter_tibble = tibble::tribble(
    ~filter               , ~n                  ,
    "missing metabolites" , nrow(has_features)  ,
    "good features"       , nrow(good_features) ,
    "good ssf"            , nrow(good_ssf)      ,
    "good range"          , nrow(good_range)    ,
    "good rank"           , nrow(good_rank)
  )
}


compare_simple_kt_pearson_plots = function(
  positive_kt,
  negative_kt,
  positive_pearson,
  negative_pearson
) {
  all_kt_other = rbind(positive_kt, negative_kt)
  all_pearson = rbind(positive_pearson, negative_pearson)
  all_pearson[is.na(all_pearson$cor), "cor"] = 0

  compare_df = data.frame(ici_kt = all_kt_other$cor, pearson = all_pearson$cor)
  set.seed(1234)
  rand_rows = sample(nrow(compare_df), 10000)
  compare_df2 = compare_df[rand_rows, ]
  out_plot = ggplot(compare_df2, aes(x = pearson, y = ici_kt)) +
    geom_point() +
    labs(x = "Pearson", y = "ICI-Kt")
  out_plot
}
