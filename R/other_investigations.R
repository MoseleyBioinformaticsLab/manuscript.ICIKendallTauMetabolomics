# these are investigations that aren't a part of the actual analysis
# per se, as in we don't want them as part of the overall workflow and
# counting towards the number of items handled by the `targets` workflow.

# checking ranges of data values
range_investigation = function(check_data) {
  # tar_load(check_data)
  low_range = check_data |>
    dplyr::filter(RANGE_CHECK %in% "TOO LOW")
  low_ids = low_range$ID

  all_to_check = paste0("checked_", low_ids)

  range_data = purrr::map(all_to_check, range_check, .progress = TRUE) |>
    purrr::list_rbind()
  check_range = dplyr::left_join(low_range, range_data, by = "ID")
  return(check_range)
}

range_check = function(checked_id) {
  # checked_id = all_to_check[20]
  checked_data = tar_read_raw(checked_id)

  ssf_data = checked_data$SUBJECT_METADATA

  measurement_data = checked_data$MEASUREMENTS[, ssf_data$sample_id] |>
    as.matrix()
  tibble::tibble(
    min = min(measurement_data, na.rm = TRUE),
    max = max(measurement_data, na.rm = TRUE),
    ID = checked_data$ID
  ) |>
    dplyr::mutate(range = max - min)
}


qc_pool_blank_investigation = function() {
  meta_data = tar_meta()
  keep_processed = meta_data |>
    dplyr::filter(grepl("^processed_", name), size != "s44b") |>
    dplyr::pull(name)

  filtered_samples = purrr::map(
    keep_processed,
    find_filtered,
    .progress = TRUE
  ) |>
    purrr::list_rbind()
}

find_filtered = function(in_processed) {
  # in_processed = "processed_AN005114"
  processed_data = tar_read_raw(in_processed)

  if (is.null(processed_data$SUBJECT_SAMPLE_FACTORS)) {
    return(NULL)
  }
  ssf_data = processed_data$SUBJECT_SAMPLE_FACTORS

  ssf_qc = ssf_data |>
    dplyr::filter(
      grepl(
        ".*pool.*|.*qc.*|.*blank.*",
        sample_id,
        ignore.case = TRUE
      ) |
        grepl(".*pool.*|.*qc.*|.*blank.*", item, ignore.case = TRUE)
    )
  ssf_qc$dataset = processed_data$ID
  ssf_qc
}
