## Load your packages, e.g. library(targets).
targets::tar_source(c("./packages.R", "R"))


# tar_option_set(
#   controller = crew_controller_local(workers = 2)
# )

# this is great, but we don't actually want to do this. We need something we can subset to
# things that are more useful.
#
# New plan:
#  save parsed data to sets of rds files that can be read in later for actual analysis;
#  figure out a way to have the mapping of samples to unique conditions returned, and
#    remove samples that aren't in groupings with enough samples;
#  do the ranking of metabolites vs num-missing and check correlations to find
#    possibly suspicious groupings (i.e. datasets that shouldn't be used)
#  save the type of measurement, and minimum per group, so we can do subsequent selections;

correlation_methods = c(
  "ici",
  "ici_completeness",
  "pearson_base",
  "pearson_base_nozero",
  "pearson_log1p",
  "pearson_log",
  "kt"
)

json_files = fs::dir_ls(
  "mwtab/repaired/json_2025-11-12"
)
mwtab_datasets = extract_mwtab_ids(json_files)

ancillary_path = "mwtab/ancillary/"

mwtab_targets = tar_map(
  mwtab_datasets,
  names = id,
  tar_target(dataset, file, format = "file"),
  tar_target(processed, parse_json(dataset, id, ancillary_path)),
  tar_target(
    checked,
    run_mwtab_checks_json(
      processed,
      min_n = 5,
      min_metabolites = 100,
      max_min_value = 20,
      min_ssf = 2,
      use_ssf_only = "yes"
    )
  ),
  tar_target(check_result, get_check(checked)),
  tar_target(smd, convert_mwtab_json_smd(checked)),
  tar_target(cor, run_cor_everyway_new(smd)),
  tar_target(limma, filter_outliers_do_limma(cor, smd)),
  tar_target(compare, limma_compare_significant(limma)),
  tar_target(missingness, calculate_missingness(smd)),
  tar_target(missingness_ranks, calculate_missingness_ranks(smd)),
  tar_target(
    missingness_ranks_correlation,
    calculate_missingness_rank_correlation(missingness_ranks)
  )
)

check_comb = tar_combine(
  check_data,
  mwtab_targets[[4]],
  command = bind_rows(!!!.x)
)

compare_comb = tar_combine(
  compare_data,
  mwtab_targets[[8]],
  command = bind_rows(!!!.x)
)

missingness_comb = tar_combine(
  missingness_data,
  mwtab_targets[[9]],
  command = bind_rows(!!!.x)
)

rank_cor_comb = tar_combine(
  rank_correlations,
  mwtab_targets[[11]],
  command = bind_rows(!!!.x)
)


mwtab_result_plan = tar_assign({
  compare_description_good = examine_limma_significant(compare_data, "GOOD") |>
    tar_target()
  compare_description_all = examine_limma_significant(compare_data, "ALL") |>
    tar_target()

  compare_stats_good = test_limma_significant(compare_data, "GOOD") |>
    tar_target()
  compare_stats_all = test_limma_significant(compare_data, "GOOD") |>
    tar_target()
})


## lod mapping as orders of magnitude ---------
lod_ranges = tibble::tibble(max = c(0.5, 1, 1.5), id = c("low", "med", "high"))

small_realistic_examples = tar_assign({
  # small theoretical example --------
  x = seq(1, 10) |>
    tar_target()
  y = seq(1, 10) |>
    tar_target()
  y2 = seq(10, 1) |>
    tar_target()
  where_na = create_na_indices(20) |>
    tar_target()

  all_kt = all_kendalltau(x, y, y2, where_na) |>
    tar_target()

  positive_kt = compare_positive_kt(x, y, where_na) |>
    tar_target()
  negative_kt = compare_negative_kt(x, y2, where_na) |>
    tar_target()
  positive_pearson = compare_positive_pearson(x, y, where_na) |>
    tar_target()
  negative_pearson = compare_negative_pearson(x, y2, where_na) |>
    tar_target()

  positive_kendall = compare_positive_pearson(
    x,
    y,
    where_na,
    method = "kendall"
  ) |>
    tar_target()
  negative_kendall = compare_negative_pearson(
    x,
    y2,
    where_na,
    method = "kendall"
  ) |>
    tar_target()

  # bigger, more realistic example ---------
  realistic_sample_1 = create_sample(n = 1000) |>
    tar_target()
  realistic_sample_2 = create_sample(n = 1000) |>
    tar_target()

  realistic_neg_sample = sort(realistic_sample_2, decreasing = TRUE) |>
    tar_target()
  realistic_neg_sample_2 = (-1 * realistic_sample_2) |>
    tar_target()
  realistic_na = create_random_na() |>
    tar_target()

  realistic_positive_kt = compare_positive_kt(
    realistic_sample_1,
    realistic_sample_2,
    realistic_na
  ) |>
    tar_target()
  realistic_negative_kt = compare_negative_kt(
    realistic_sample_1,
    realistic_neg_sample,
    realistic_na
  ) |>
    tar_target()
  realistic_negative_kt_2 = compare_negative_kt(
    realistic_sample_1,
    realistic_neg_sample_2,
    realistic_na
  ) |>
    tar_target()

  realistic_positive_pearson = compare_positive_pearson(
    realistic_sample_1,
    realistic_sample_2,
    realistic_na
  ) |>
    tar_target()
  realistic_negative_pearson = compare_negative_pearson(
    realistic_sample_1,
    realistic_neg_sample,
    realistic_na
  ) |>
    tar_target()
  realistic_negative_pearson_2 = compare_negative_pearson(
    realistic_sample_1,
    realistic_neg_sample_2,
    realistic_na
  ) |>
    tar_target()

  realistic_positive_kendall = compare_positive_kendall(
    realistic_sample_1,
    realistic_sample_2,
    realistic_na
  ) |>
    tar_target()
  realistic_negative_kendall = compare_negative_kendall(
    realistic_sample_1,
    realistic_neg_sample,
    realistic_na
  ) |>
    tar_target()

  realistic_reference_plot = compare_realistic_to_reference(
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
  ) |>
    tar_target()
  realistic_comparison_plot = compare_realistic_to_each(
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
  ) |>
    tar_target()

  left_censored_samples = create_lc_samples() |>
    tar_target()
  left_censored_cor = left_censor_correlate(left_censored_samples) |>
    tar_target()
  random_censored_cor = random_censor_correlate(left_censored_samples) |>
    tar_target()
  logtransform_censored_cor = lt_left_censor_correlate(left_censored_samples) |>
    tar_target()

  subsample = sample(1000, 50) |>
    tar_target()
  left_sample_cor = left_censor_correlate(left_censored_samples[subsample, ]) |>
    tar_target()
  random_sample_cor = random_censor_correlate(
    left_censored_samples[subsample, ],
    n_na = seq(0, 12, 2)
  ) |>
    tar_target()
  logtransform_sample_cor = lt_left_censor_correlate(left_censored_samples[
    subsample,
  ]) |>
    tar_target()

  censored_value_plots = plot_censored_data(
    left_censored_cor,
    random_censored_cor,
    logtransform_censored_cor
  ) |>
    tar_target()
  censored_compare_plots = compare_censored_data(
    left_censored_cor,
    random_censored_cor,
    logtransform_censored_cor
  ) |>
    tar_target()
})

vl_plan = tar_assign({
  # variable lod data creation --------
  lod_ranges_tar = lod_ranges |>
    tar_target()
  check_lod_levels = seq(0.1, 3, by = 0.1) |>
    tar_target()
  lod_vars = list(
    n_feature = 1000,
    n_sample = 100,
    meanlog = 1,
    sdlog = 0.5,
    sd = 0.2
  ) |>
    tar_target()
  var_lod_samples = create_large_replicate_samples(lod_vars) |>
    tar_target()

  ## verify how the changes in missing some order of magnitude introduces missing values
  vl_na_perc = check_lod_na_perc(var_lod_samples, check_lod_levels) |>
    tar_target()

  vl_na_perc_graph = create_na_perc_graph(vl_na_perc, lod_ranges_tar) |>
    tar_target()
  vl_diff_graph = create_lod_diff_graph(vl_cor_diff_all) |>
    tar_target()
  vl_icikt_kt_graph = create_icikt_ktimpute_graph(vl_cor_diff_all) |>
    tar_target()
})

vl_lod_map = tar_map(
  lod_ranges,
  names = id,
  tar_target(vl_samples, create_variable_lod_samples(var_lod_samples, max, id)),
  tar_target(vl_cor, calculate_variable_correlations(vl_samples)),
  tar_target(vl_cor_diff, calculate_var_lod_correlation_diffs(vl_cor))
)
# tar_target(vl_diff_summary,
#            calculate_cor_diff_summaries(vl_cor_diff)))
#
vl_cor_diff_combine_map = tar_combine(
  vl_cor_diff_all,
  vl_lod_map[[3]],
  command = bind_rows(!!!.x)
)


figures_tables_plan = tar_assign({
  # ---- missingness figures
  missingness_summary = missingness_data |>
    dplyr::mutate(padjust = p.adjust(p.value, method = "BH")) |>
    tar_target()
  missingness_percentage_sina_plot = create_missing_percentage_sina_plot() |>
    tar_target()

  variable_dynamic_range_image = create_variable_dynamic_range_image(
    vl_diff_graph
  ) |>
    tar_target()
  censored_value_plot_image = create_censored_value_images(
    censored_value_plots
  ) |>
    tar_target()
  rank_missingness_image = create_figure1_image(
    missingness_summary,
    missingness_tests_AN001074,
    rank_correlations
  ) |>
    tar_target()
})

docs_plan = tar_assign({
  ici_kt_manuscript = tar_render("docs/ici_kt_manuscript.Rmd")
  supp_materials = tar_render("docs/supplemental_materials.Rmd")
})

list(
  mwtab_targets,
  check_comb,
  compare_comb,
  mwtab_result_plan,
  missingness_comb,
  rank_cor_comb,
  small_realistic_examples,
  vl_plan,
  vl_lod_map,
  vl_cor_diff_combine_map,
  figures_tables_plan
)
