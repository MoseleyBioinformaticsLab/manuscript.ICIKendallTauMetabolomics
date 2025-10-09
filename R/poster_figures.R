create_example_missingness = function(in_cor) {
  # in_cor = tar_read("metabolomics_cor_NSCLC")
  use_counts = assays(in_cor$data)$normalized
  use_info = colData(in_cor$data) |> tibble::as_tibble()
  rank_counts = ICIKendallTau::rank_order_data(
    use_counts,
    sample_classes = use_info$factors
  )

  rank_images = purrr::map(rank_counts, \(in_counts) {
    use_values = in_counts$ordered |>
      as.data.frame()
    visdat::vis_miss(use_values, show_perc = FALSE) +
      scale_x_discrete(labels = NULL, breaks = NULL) +
      labs(x = NULL) +
      theme(legend.position = "none")
  })

  basic_ragg("docs/poster/working/missing_example.png")
  print(rank_images[[1]])
  dev.off()
}


rank_missingness = function(in_cor) {
  # in_cor = tar_read("metabolomics_cor_AN001074")
  # in_cor = tar_read("metabolomics_cor_NSCLC")
  use_ranks = in_cor$ranks
  median_image = use_ranks |>
    ggplot(aes(x = n_na, y = median_rank)) +
    geom_point(size = 2) +
    facet_wrap(~factors) +
    labs(x = "Number of Samples Missing In", y = "Median Rank") +
    theme(strip.text = element_text(size = 10))

  use_ranks_min = use_ranks |>
    dplyr::summarise(min_rank = min(median_rank), .by = c(factors, n_na))

  min_image = use_ranks_min |>
    ggplot(aes(x = n_na, y = min_rank)) +
    geom_point(size = 2) +
    facet_wrap(~factors) +
    labs(x = "Number of Samples Missing In", y = "Min(Median Rank)") +
    theme(strip.text = element_text(size = 10))

  basic_ragg("docs/poster/images/median_ranks.png")
  print(median_image)
  dev.off()

  basic_ragg("docs/poster/images/min_ranks.png")
  print(min_image)
  dev.off()

  return(c("median_ranks.png", "min_ranks.png"))
}

create_percent_missingness = function(missingness_summary) {
  # tar_load(missingness_summary)
  out_plot = missingness_summary |>
    ggplot(aes(x = perc_miss, y = type)) +
    geom_sina(size = 2) +
    labs(x = "% Missing Value", y = "Method")

  basic_ragg("docs/poster/images/missingness_percentage.png")
  print(out_plot)
  dev.off()
}

create_rank_correlation_image = function(rank_correlations) {
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
      )
    ) |>
    ggplot(aes(x = correlation, y = rank_type)) +
    geom_sina(size = 2) +
    facet_wrap(~type, ncol = 1) +
    labs(x = "Kendall-tau(Rank vs N-Missing)", y = "Rank Summary")

  basic_ragg("docs/poster/images/rank_correlation_distribution.png")
  print(out_plot)
  dev.off()
  return(invisible(NULL))
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
