check_intersect_analyses = function() {
  enrichment_datasets_dir = "/big_data/data/metabolite_predictions/PEAManuscriptData/PathwayEnrichmentAnalysis/data/enrichment/annotations"

  check_res = tar_read(check_data)
  check_good = (check_res$FEATURE_CHECK %in% "GOOD") &
    (check_res$SSF_CHECK %in% "GOOD SSF") &
    (check_res$RANK_CHECK %in% c("GOOD", "SIGN DIFFERENCE")) &
    (check_res$RANGE_CHECK %in% "GOOD")

  check_res_good = check_res[check_good, ]

  enrichment_datasets_list = tibble::tibble(
    dir = fs::dir_ls(enrichment_datasets_dir)
  ) |>
    dplyr::mutate(ID = fs::path_file(dir))
}

check_n_reps = function() {
  enrichment_datasets_dir = "/big_data/data/metabolite_predictions/PEAManuscriptData/PathwayEnrichmentAnalysis/data/enrichment/annotations"

  meta_data = tar_meta(targets_only = TRUE)
  smd_good = meta_data |>
    dplyr::filter(grepl("^smd", name) & !(size == "s44b"))

  enrichment_datasets_list = tibble::tibble(
    dir = fs::dir_ls(enrichment_datasets_dir)
  ) |>
    dplyr::mutate(ID = fs::path_file(dir))

  smd_good$ID = gsub("smd_", "", smd_good$name)
  smd_intersect = dplyr::inner_join(
    smd_good,
    enrichment_datasets_list,
    by = "ID"
  )

  data_dir = "/big_data/data/rmflight_icikt_metabolomics/_targets/objects"
  n_reps_features = purrr::map(smd_intersect$name, \(in_file) {
    tmp_data = readRDS(fs::path(data_dir, in_file))
    tibble::tibble(n_samples = ncol(tmp_data), n_features = nrow(tmp_data))
  }) |>
    purrr::list_rbind()
  smd_intersect = dplyr::bind_cols(smd_intersect, n_reps_features)
}

check_mapping_ids = function(smd, grouped_annotations) {
  smd = tar_read(smd_AN002783)
  grouped_annotations = tar_read(predicted_annotations_grouped)

  match_id = metadata(smd)$CHECK$ID
  match_annotation = grouped_annotations |>
    dplyr::filter(id %in% match_id)

  feature_data = rowData(smd) |>
    tibble::as_tibble() |>
    dplyr::mutate(is_smd = TRUE)

  annotation_features = match_annotation |>
    dplyr::select(features) |>
    dplyr::distinct() |>
    dplyr::mutate(is_annotation = TRUE)
  match_features = dplyr::full_join(
    feature_data,
    annotation_features,
    by = c("metabolite" = "features")
  )
  match_features |>
    dplyr::select(metabolite, is_smd, is_annotation) |>
    View()
}
