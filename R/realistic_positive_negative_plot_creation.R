compare_realistic_to_reference = function(
  realistic_sample_1,
  realistic_sample_2,
  realistic_neg_sample_2,
  realistic_na,
  realistic_positive_pearson,
  realistic_positive_kendall,
  realistic_positive_kt,
  realistic_negative_pearson_2,
  realistic_negative_kendall,
  realistic_negative_kt_2
) {
  # tar_load(realistic_sample_1)
  # tar_load(realistic_sample_2)
  # tar_load(realistic_neg_sample_2)
  n_na = purrr::map_int(realistic_na, length)
  ref_pearson = cor(realistic_sample_1, realistic_sample_2)
  ref_kendall = ici_kt(realistic_sample_1, realistic_sample_2, "global")[1]

  ref_pearson_neg = cor(realistic_sample_1, realistic_neg_sample_2)
  ref_kendall_neg = ici_kt(
    realistic_sample_1,
    realistic_neg_sample_2,
    "global"
  )[1]

  # tar_load(realistic_na)
  # tar_load(realistic_positive_pearson)
  # tar_load(realistic_positive_kendall)
  # tar_load(realistic_positive_kt)
  # tar_load(realistic_negative_pearson_2)
  # tar_load(realistic_negative_kendall)
  # tar_load(realistic_negative_kt_2)

  realistic_positive_pearson = realistic_positive_pearson %>%
    dplyr::mutate(
      type = "Pearson",
      dir = "positive",
      diff = ref_pearson - cor,
      n_na = (x_na + y_na) / 2
    )

  realistic_negative_pearson_2 = realistic_negative_pearson_2 %>%
    dplyr::mutate(
      type = "Pearson",
      dir = "negative",
      diff = ref_pearson_neg - cor,
      n_na = (x_na + y_na) / 2
    )

  realistic_positive_kendall = realistic_positive_kendall %>%
    dplyr::mutate(
      type = "Kendall",
      dir = "positive",
      diff = ref_kendall - cor,
      n_na = (x_na + y_na) / 2
    )

  realistic_negative_kendall = realistic_negative_kendall %>%
    dplyr::mutate(
      type = "Kendall",
      dir = "negative",
      diff = ref_kendall_neg - cor,
      n_na = (x_na + y_na) / 2
    )

  realistic_positive_kt = realistic_positive_kt %>%
    dplyr::mutate(
      type = "ICI-Kt",
      dir = "positive",
      diff = ref_kendall - cor,
      n_na = (x_na + y_na) / 2
    )

  realistic_negative_kt_2 = realistic_negative_kt_2 %>%
    dplyr::mutate(
      type = "ICI-Kt",
      dir = "negative",
      diff = ref_kendall_neg - cor,
      n_na = (x_na + y_na) / 2
    )

  positive_df = rbind(
    realistic_positive_pearson,
    realistic_positive_kendall,
    realistic_positive_kt
  )

  negative_df = rbind(
    realistic_negative_pearson_2,
    realistic_negative_kendall,
    realistic_negative_kt_2
  )

  use_rows = rep(FALSE, nrow(negative_df))
  n_subset = 10000
  use_rows[sample(length(use_rows), n_subset)] = TRUE
  negative_df = negative_df %>%
    dplyr::mutate(
      perc_missing = n_na / max(n_na) * 100,
      type = factor(
        type,
        levels = c("Kendall", "Pearson", "ICI-Kt"),
        ordered = TRUE
      )
    ) %>%
    dplyr::filter(use_rows)
  positive_df = positive_df %>%
    dplyr::mutate(
      perc_missing = n_na / max(n_na) * 100,
      type = factor(
        type,
        levels = c("Kendall", "Pearson", "ICI-Kt"),
        ordered = TRUE
      )
    ) %>%
    dplyr::filter(use_rows)
  pos_diff = ggplot(positive_df, aes(x = n_na, y = diff, color = n_na)) +
    geom_point() +
    facet_grid(type ~ ., scales = "free") +
    labs(
      subtitle = "Positive Correlation",
      x = "Number Missing",
      y = "Difference from Reference of 1",
      color = "# NA"
    ) +
    scale_color_viridis_c() +
    scale_y_continuous(labels = label_scientific(digits = 2)) +
    theme(plot.subtitle = element_text(size = 13))
  neg_diff = ggplot(negative_df, aes(x = n_na, y = diff, color = n_na)) +
    geom_point(show.legend = FALSE) +
    scale_color_viridis_c() +
    facet_grid(type ~ ., scales = "free") +
    labs(
      subtitle = "Negative Correlation",
      x = "Number Missing",
      y = "Difference from Reference of -1",
      color = "# NA"
    ) +
    scale_y_continuous(labels = label_scientific(digits = 2)) +
    theme(plot.subtitle = element_text(size = 13))

  list(positive = pos_diff, negative = neg_diff, n_subset = n_subset)
}

compare_realistic_to_each = function(
  realistic_sample_1,
  realistic_sample_2,
  realistic_neg_sample_2,
  realistic_na,
  realistic_positive_pearson,
  realistic_positive_kendall,
  realistic_positive_kt,
  realistic_negative_pearson_2,
  realistic_negative_kendall,
  realistic_negative_kt_2
) {
  # tar_load(realistic_sample_1)
  # tar_load(realistic_sample_2)
  # tar_load(realistic_neg_sample_2)
  # tar_load(realistic_na)
  # tar_load(realistic_positive_pearson)
  # tar_load(realistic_positive_kendall)
  # tar_load(realistic_positive_kt)
  # tar_load(realistic_negative_pearson_2)
  # tar_load(realistic_negative_kendall)
  # tar_load(realistic_negative_kt_2)

  # ref_pearson = cor(realistic_sample_1, realistic_sample_2)
  # ref_kendall = ici_kt(realistic_sample_1, realistic_sample_2, "global")[1]
  #
  # ref_pearson_neg = cor(realistic_sample_1, realistic_neg_sample)
  # ref_kendall_neg = ici_kt(realistic_sample_1, realistic_neg_sample, "global")[1]

  # set.seed(0134)
  n_na = purrr::map_int(realistic_na, length)

  rp_pearson_wide = realistic_positive_pearson |>
    dplyr::transmute(
      pearson = cor,
      id = paste0(i_na, ".", x_na, ".", y_na),
      x_na = x_na,
      y_na = y_na,
      n_na = (x_na + y_na) / 2
    )

  rp_kendall_wide = realistic_positive_kendall |>
    dplyr::transmute(kendall = cor, id = paste0(i_na, ".", x_na, ".", y_na))

  rp_icikt_wide = realistic_positive_kt |>
    dplyr::transmute(icikt = cor, id = paste0(i_na, ".", x_na, ".", y_na))

  compare_positive = dplyr::left_join(
    rp_pearson_wide,
    rp_kendall_wide,
    by = "id"
  )
  compare_positive = dplyr::left_join(
    compare_positive,
    rp_icikt_wide,
    by = "id"
  )

  # use_rows = rep(TRUE, nrow(compare_positive))
  use_rows = rep(FALSE, nrow(compare_positive))
  n_subset = 10000
  keep_rows = sample(length(use_rows), n_subset)
  use_rows[keep_rows] = TRUE

  compare_positive_subset = compare_positive[use_rows, ]
  compare_positive_subset = compare_positive_subset |>
    dplyr::rowwise() |>
    dplyr::mutate(max_na = max(c(x_na, y_na))) |>
    dplyr::ungroup()

  rp_pearson_negative_wide = realistic_negative_pearson_2 %>%
    dplyr::transmute(
      pearson = cor,
      x_na = x_na,
      y_na = y_na,
      id = paste0(i_na, ".", x_na, ".", y_na),
      n_na = (x_na + y_na) / 2
    )
  rp_kendall_negative_wide = realistic_negative_kendall |>
    dplyr::transmute(kendall = cor, id = paste0(i_na, ".", x_na, ".", y_na))
  rp_icikt_negative_wide = realistic_negative_kt_2 |>
    dplyr::transmute(icikt = cor, id = paste0(i_na, ".", x_na, ".", y_na))

  compare_negative = dplyr::left_join(
    rp_pearson_negative_wide,
    rp_kendall_negative_wide,
    by = "id"
  )
  compare_negative = dplyr::left_join(
    compare_negative,
    rp_icikt_negative_wide,
    by = "id"
  )

  compare_negative_subset = compare_negative[use_rows, ]
  compare_negative_subset = compare_negative_subset |>
    dplyr::rowwise() |>
    dplyr::mutate(max_na = max(c(x_na, y_na))) |>
    dplyr::ungroup()

  p_ici_pos = compare_positive_subset |>
    ggplot(aes(x = pearson, y = icikt, color = n_na)) +
    scale_color_viridis_c() +
    geom_point() +
    labs(
      x = "Pearson",
      y = "ICI-Kt",
      subtitle = "Positive Correlation: 1",
      color = "# NA"
    ) +
    theme(legend.position = c(0.2, 0.8))
  k_ici_pos = compare_positive_subset |>
    ggplot(aes(x = kendall, y = icikt, color = n_na)) +
    scale_color_viridis_c() +
    geom_point(show.legend = FALSE) +
    labs(x = "Kendall", y = "ICI-Kt")

  p_ici_neg = compare_negative_subset |>
    ggplot(aes(x = pearson, y = icikt, color = n_na)) +
    scale_color_viridis_c() +
    geom_point(show.legend = FALSE) +
    labs(x = "Pearson", y = "ICI-Kt", subtitle = "Negative Correlation: -1")
  k_ici_neg = compare_negative_subset |>
    ggplot(aes(x = kendall, y = icikt, color = n_na)) +
    geom_point(show.legend = FALSE) +
    scale_color_viridis_c() +
    labs(x = "Kendall", y = "ICI-Kt")

  list(
    positive = list(pearson = p_ici_pos, kendall = k_ici_pos),
    negative = list(pearson = p_ici_neg, kendall = k_ici_neg),
    n_subset = n_subset
  )
}

plot_censored_data = function(
  left_censored_cor,
  random_censored_cor,
  logtransform_censored_cor
) {
  # tar_load(left_censored_cor)
  # tar_load(random_censored_cor)
  # tar_load(logtransform_censored_cor)

  left_censored_cor = left_censored_cor %>%
    dplyr::mutate(
      which2 = dplyr::case_when(
        which %in% "ici" ~ "ICI-Kt",
        which %in% "kendall" ~ "Kendall",
        which %in% "kendall_0" ~ "Kendall-0",
        which %in% "pearson" ~ "Pearson",
        which %in% "pearson_0" ~ "Pearson-0"
      )
    )

  left_ref_y = left_censored_cor |>
    dplyr::group_by(which) |>
    dplyr::summarise(ref_y = max(cor))
  left_ref_cor = left_censored_cor |>
    dplyr::filter(n_na == 0) |>
    dplyr::group_by(which) |>
    dplyr::summarise(ref_cor = cor[1]) |>
    dplyr::mutate(
      ref_str = format(ref_cor, digits = 3),
      ref_x = dplyr::case_match(
        which,
        "ici" ~ 25,
        "kendall" ~ 100,
        "kendall_0" ~ 25,
        "pearson" ~ 100,
        "pearson_0" ~ 100
      )
    )
  left_ref_cor = dplyr::left_join(left_ref_cor, left_ref_y, by = "which")
  left_ref_cor = dplyr::left_join(
    left_ref_cor,
    left_censored_cor |>
      dplyr::select(which, which2) |>
      dplyr::distinct(),
    by = "which"
  )

  random_censored_cor = random_censored_cor %>%
    dplyr::mutate(censoring = "random") %>%
    dplyr::filter(n_na > 0) %>%
    dplyr::mutate(
      which2 = dplyr::case_when(
        which %in% "ici" ~ "ICI-Kt",
        which %in% "kendall" ~ "Kendall",
        which %in% "kendall_0" ~ "Kendall-0",
        which %in% "pearson" ~ "Pearson",
        which %in% "pearson_0" ~ "Pearson-0"
      )
    )

  random_ref_cor = left_ref_cor |>
    dplyr::arrange(which) |>
    dplyr::mutate(
      ref_x = c(150, 25, 150, 25, 150),
      ref_y = c(0.8, 0.880, 0.8, 0.99, 0.9)
    )

  logtransform_censored_cor = logtransform_censored_cor %>%
    dplyr::mutate(censoring = "log") %>%
    dplyr::filter(n_na > 0) %>%
    dplyr::mutate(
      which2 = dplyr::case_when(
        which %in% "ici" ~ "ICI-Kt",
        which %in% "kendall" ~ "Kendall",
        which %in% "kendall_0" ~ "Kendall-0",
        which %in% "pearson" ~ "Pearson",
        which %in% "pearson_0" ~ "Pearson-0"
      )
    )

  logtransform_ref_y = logtransform_censored_cor |>
    dplyr::group_by(which) |>
    dplyr::summarise(ref_y = max(cor))
  logtransform_ref_cor = logtransform_censored_cor |>
    dplyr::filter(n_na == 2) |>
    dplyr::group_by(which) |>
    dplyr::summarise(ref_cor = cor[1]) |>
    dplyr::mutate(
      ref_str = format(ref_cor, digits = 3),
      ref_x = dplyr::case_match(
        which,
        "ici" ~ 25,
        "kendall" ~ 100,
        "kendall_0" ~ 25,
        "pearson" ~ 25,
        "pearson_0" ~ 100
      )
    )
  logtransform_ref_cor = dplyr::left_join(
    logtransform_ref_cor,
    logtransform_ref_y,
    by = "which"
  )
  logtransform_ref_cor = dplyr::left_join(
    logtransform_ref_cor,
    logtransform_censored_cor |>
      dplyr::select(which, which2) |>
      dplyr::distinct(),
    by = "which"
  )

  random_plot = ggplot(random_censored_cor, aes(x = n_na, y = cor)) +
    geom_sina(aes(group = n_na)) +
    geom_text(
      data = random_ref_cor,
      aes(x = ref_x, y = ref_y, label = ref_str),
      hjust = "left"
    ) +
    facet_wrap(~which2, nrow = 1, scales = "free_y") +
    labs(x = "Number Missing", y = "Correlation") +
    cowplot::panel_border()

  rp_build = ggplot_build(random_plot)
  rp_breaks = rp_build$layout$panel_params[[1]]$x$get_breaks()
  rp_lim = rp_build$layout$panel_params[[1]]$x$get_limits()

  left_plot = ggplot(left_censored_cor, aes(x = n_na, y = cor)) +
    geom_point() +
    geom_text(
      data = left_ref_cor,
      aes(x = ref_x, y = ref_y, label = ref_str),
      hjust = "left"
    ) +
    facet_wrap(~which2, nrow = 1, scales = "free_y") +
    labs(x = "Number Missing", y = "Correlation") +
    scale_x_continuous(breaks = c(0, 100, 200)) +
    cowplot::panel_border()

  log_plot = ggplot(logtransform_censored_cor, aes(x = n_na, y = cor)) +
    geom_point() +
    geom_text(
      data = logtransform_ref_cor,
      aes(x = ref_x, y = ref_y, label = ref_str),
      hjust = "left"
    ) +
    facet_wrap(~which2, nrow = 1, scales = "free_y") +
    labs(x = "Number Missing", y = "Correlation") +
    scale_x_continuous(breaks = c(0, 100, 200)) +
    cowplot::panel_border()
  list(left = left_plot, log = log_plot, random = random_plot)
}

compare_censored_data = function(
  left_censored_cor,
  random_censored_cor,
  logtransform_censored_cor
) {
  # tar_load(left_censored_cor)
  # tar_load(random_censored_cor)
  # tar_load(logtransform_censored_cor)

  left_plots = comparison_individual_plot_maker(left_censored_cor, tag = "A")
  log_plots = comparison_individual_plot_maker(
    logtransform_censored_cor,
    tag = "B"
  )
  random_censored_cor = random_censored_cor |>
    dplyr::filter(n_na > 0) |>
    dplyr::mutate(id = paste0(n_na, "_", cutoff, "_", rep))
  random_plots = comparison_individual_plot_maker(
    random_censored_cor,
    tag = "C"
  )

  list(left = left_plots, log = log_plots, random = random_plots)
}

comparison_individual_plot_maker = function(censored_cor, tag = "") {
  censored_cor = censored_cor |>
    dplyr::mutate(
      which2 = dplyr::case_when(
        which %in% "ici" ~ "ICI-Kt",
        which %in% "kendall" ~ "Kendall",
        which %in% "kendall_0" ~ "Kendall-0",
        which %in% "pearson" ~ "Pearson",
        which %in% "pearson_0" ~ "Pearson-0"
      )
    )

  if (is.null(censored_cor$id)) {
    censored_cor = censored_cor |>
      dplyr::mutate(id = paste0(n_na, "_", cutoff))
  }
  censored_wide = censored_cor |>
    tidyr::pivot_wider(id_cols = id, names_from = "which2", values_from = "cor")
  censored_pearson = censored_wide |>
    ggplot(aes(x = `Pearson`, y = `ICI-Kt`)) +
    geom_point() +
    labs(tag = tag)
  censored_pearson_0 = censored_wide |>
    ggplot(aes(x = `Pearson-0`, `ICI-Kt`)) +
    geom_point()
  censored_kendall = censored_wide |>
    ggplot(aes(x = `Kendall`, y = `ICI-Kt`)) +
    geom_point()
  censored_kendall_0 = censored_wide |>
    ggplot(aes(x = `Kendall-0`, y = `ICI-Kt`)) +
    geom_point()
  list(
    pearson = censored_pearson,
    pearson_0 = censored_pearson_0,
    kendall = censored_kendall,
    kendall_0 = censored_kendall_0
  )
}
