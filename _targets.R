## Load your packages, e.g. library(targets).
targets::tar_source(c("./packages.R", "R"))


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


correlation_samples = readRDS("data/smd/mwtab_smd.rds")


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

metabolomics_map = tar_map(
  correlation_samples,
  names = id,
  tar_target(
    metabolomics_cor,
    run_cor_everyway_new(id, smd_file, type, value_check)
  ),
  tar_target(
    outliers,
    calculate_outlier_effects(metabolomics_cor)
  )
)

list(
  small_realistic_examples,
  vl_plan,
  vl_lod_map,
  vl_cor_diff_combine_map,
  metabolomics_map
)
