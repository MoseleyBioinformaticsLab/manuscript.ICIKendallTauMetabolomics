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

create_parallel_performance_figure = function(
  test_results = res_12,
  log_dir = "python"
) {
  log_files = fs::dir_ls(log_dir, regexp = "r_log")
  res_files = fs::dir_ls(log_dir, regexp = "r_res")

  res_data = purrr::map(res_files, \(in_file) {
    # in_file = res_files[1]
    readLines(in_file) |>
      stringr::str_replace_all('"', '') |>
      tibble::as_tibble_row(.name_repair = "unique")
  }) |>
    purrr::list_rbind()

  names(res_data) = c("cores", "time", "log")
  res_data$cores = as.numeric(res_data$cores)
  res_data$time = as.numeric(res_data$time)

  time_plot = res_data |>
    ggplot(aes(x = cores, y = time)) +
    geom_line() +
    geom_point() +
    theme(legend.position = "inside", legend.position.inside = c(0.2, 0.8)) +
    labs(x = "# Cores", y = "Time(s)")

  log_data = purrr::map(log_files, \(in_file) {
    # in_file = log_files[1]
    log_lines = readLines(in_file)
    used_mem = stringr::str_extract(log_lines, "Active: [0-9]+") |>
      stringr::str_replace("Active: ", "") |>
      as.numeric()

    used_mem = used_mem[!is.na(used_mem)]
    tibble::tibble(log = in_file, memory = (max(used_mem) - used_mem[1]) / 1032)
  }) |>
    purrr::list_rbind()

  res_memory = dplyr::left_join(res_data, log_data, by = "log")

  memory_plot = res_memory |>
    ggplot(aes(x = cores, y = memory)) +
    geom_line() +
    geom_point() +
    labs(x = "# Cores", y = "Memory (MiB)")

  time_plot | memory_plot
}
