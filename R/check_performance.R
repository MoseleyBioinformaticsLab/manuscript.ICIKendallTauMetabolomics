create_performance_figure = function(single_core_perf) {
  # tar_load(single_core_perf)
  pearson = single_core_perf %>%
    dplyr::filter(method %in% "pearson")

  p_comp = ggplot(pearson, aes(x = n, y = median)) +
    geom_line() +
    geom_smooth(method = "lm", formula = y ~ (x)) +
    labs(
      subtitle = 'stats::cor(method = "pearson"), O(n)',
      x = "Number of Features",
      y = "Time (ns)"
    )

  kendall = single_core_perf %>%
    dplyr::filter(method %in% "kendall")

  k_comp = ggplot(kendall, aes(x = n, y = median)) +
    geom_line() +
    geom_smooth(method = "lm", formula = y ~ I(x^2)) +
    labs(
      subtitle = 'stats::cor(method = "kendall"), O(n^2)',
      x = "Number of Features",
      y = "Time (ns)"
    )

  ici = single_core_perf %>%
    dplyr::filter(method %in% "icikt")

  ici_comp = ggplot(ici, aes(x = n, y = median)) +
    geom_line() +
    geom_smooth(method = "lm", formula = y ~ x * log(x)) +
    labs(
      subtitle = "ICIKendallTau::ici_kt(), O(nlog(n))",
      x = "Number of Features",
      y = "Time (ns)"
    )

  p_comp / ici_comp / k_comp
}
