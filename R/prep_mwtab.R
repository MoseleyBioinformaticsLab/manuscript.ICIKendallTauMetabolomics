extract_mwtab_ids = function(mwtab_files) {
  split_files = stringr::str_split(fs::path_file(mwtab_files), "%2F")
  mwtab_df = tibble::tibble(
    id = purrr::map_chr(split_files, \(in_file) {
      in_file[7]
    }),
    file = mwtab_files
  )
  mwtab_df
}


parse_mwtab = function(mwtab_file, id, ancillary_path) {
  #mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000023%2Fmwtab%2Ftxt.txt"  # MS file
  # mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000693%2Fmwtab%2Ftxt.txt"  # NMR file

  # mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000001%2Fmwtab%2Ftxt.txt"

  #mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000004%2Fmwtab%2Ftxt.txt"

  # mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000038%2Fmwtab%2Ftxt.txt"

  # mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000555%2Fmwtab%2Ftxt.txt"

  # mwtab_file = tar_read(mwtab_AN000172)
  # mwtab_file = tar_read(mwtab_AN003177)
  # mwtab_file = tar_read(mwtab_AN002445)
  # mwtab_file = tar_read(mwtab_AN003894)
  # mwtab_file = tar_read(mwtab_AN001238)
  mwtab_lines = readLines(mwtab_file)
  block_lines = which(stringr::str_detect(
    mwtab_lines,
    "^\\#[[:upper:]]|EXTENDED.*START"
  ))
  end_line = which(stringr::str_detect(mwtab_lines, "^\\#END"))
  block_lines = block_lines[!(block_lines %in% end_line)]

  block_data = vector("list", length(block_lines))
  for (iblock in seq_len(length(block_lines))) {
    # message(iblock)
    # iblock = 1
    if (iblock != length(block_lines)) {
      block_range = seq(block_lines[iblock] + 1, block_lines[iblock + 1] - 1)
    } else {
      block_range = seq(block_lines[iblock] + 1, length(mwtab_lines) - 1)
    }
    if (iblock == 1) {
      block_id = "MWINFO"
      block_range = c(1, block_range)
    } else {
      block_id = stringr::str_replace(
        mwtab_lines[block_lines[iblock]],
        "\\#",
        ""
      ) |>
        stringr::str_replace("\\ .*", "") |>
        stringr::str_replace("\\:", "")
    }
    names(block_data)[iblock] = block_id
    block_data[[iblock]] = mwtab_lines[block_range]
  }

  parsed_data = purrr::imap(block_data, \(block, block_id) {
    switch(
      block_id,
      MWINFO = parse_tabs(block),
      PROJECT = parse_tabs(block),
      STUDY = parse_tabs(block),
      SUBJECT = parse_tabs(block),
      SUBJECT_SAMPLE_FACTORS = parse_factors_mwtab(block),
      COLLECTION = parse_tabs(block),
      TREATMENT = parse_tabs(block),
      SAMPLEPREP = parse_tabs(block),
      CHROMATOGRAPHY = parse_tabs(block),
      ANALYSIS = parse_tabs(block),
      MS = parse_tabs(block),
      MS_METABOLITE_DATA = parse_ms_data_mwtab(block),
      METABOLITES = parse_metabolites_mwtab(block),
      NMR = parse_tabs(block),
      NMR_BINNED_DATA = parse_nmr_data_mwtab(block),
      block
    )
  })

  parsed_data$ID = id
  if (is.null(parsed_data$SUBJECT_SAMPLE_FACTORS)) {
    # browser()
    return(parsed_data)
  }
  n_reps_factors = count_factors_replicates_mwtab(
    parsed_data$SUBJECT_SAMPLE_FACTORS
  )
  parsed_data$SUBJECT_SAMPLE_FACTORS = n_reps_factors$sample_factors
  parsed_data$N_REPLICATES = n_reps_factors$n

  if (!is.null(parsed_data$MS_METABOLITE_DATA)) {
    parsed_data$MEASUREMENTS = parsed_data$MS_METABOLITE_DATA
    parsed_data$MEASUREMENT_TYPE = "MS"
  } else if (!is.null(parsed_data$NMR_BINNED_DATA)) {
    parsed_data$MEASUREMENTS = parsed_data$NMR_BINNED_DATA
    parsed_data$MEASUREMENT_TYPE = "NMR"
  } else {
    parsed_data = check_and_parse_ancillary_mwtab(parsed_data)
  }

  parsed_data
}

check_and_parse_ancillary_mwtab = function(
  parsed_data,
  ancillary_path = "/big_data/data/rmflight_icikt_poster_iecm/data/ancillary/"
) {
  # parsed_data = tar_read(processed_AN000111)
  # ancillary_path = "/big_data/data/rmflight_icikt_poster_iecm/data/ancillary/"
  if (!is.null(parsed_data$MS)) {
    results_file = parsed_data$MS |>
      dplyr::filter(grepl("MS_RESULTS_FILE", field)) |>
      dplyr::pull(value)
    data_type = "MS"
  } else if (!is.null(parsed_data$NMR)) {
    results_file = parsed_data$NMR |>
      dplyr::filter(grepl("NMR_RESULTS_FILE", field)) |>
      dplyr::pull(value)
    data_type = "NMR"
  }
  just_file = stringr::str_split(results_file, "\t")[[1]][1]

  ancillary_files = fs::dir_ls(ancillary_path)
  ancillary_just_file = fs::path_file(ancillary_files)
  names(ancillary_files) = ancillary_just_file

  match_file = ancillary_files[just_file]

  # can't find a match, just return early
  if (is.na(match_file)) {
    return(parsed_data)
  }

  ancillary_data = readr::read_delim(
    match_file,
    delim = "\t",
    show_col_types = FALSE
  ) |>
    janitor::clean_names()

  check_character_reread = function(ancillary_data, match_file) {
    character_cols = purrr::map_lgl(ancillary_data, \(x) {
      inherits(x, "character")
    })
    numeric_cols = !character_cols

    if (sum(character_cols) == ncol(ancillary_data)) {
      ancillary_tmp = readr::read_delim(
        match_file,
        skip = 1,
        col_names = FALSE,
        show_col_types = FALSE
      )
      character_cols = purrr::map_lgl(ancillary_tmp, \(x) {
        inherits(x, "character")
      })
      if (sum(character_cols) == ncol(ancillary_tmp)) {
        ancillary_tmp = readr::read_delim(
          match_file,
          skip = 2,
          col_names = FALSE,
          show_col_types = FALSE
        )
      }
    }
  }
  # not everything can be a character column, there should
  # be some numeric columns

  # we assume the first character column is an ID.
  # if there are no character columns, we just make our own.
  if (length(character_cols) == 0) {
    metabolite_data = tibble::tibble(
      feature_id = janitor::make_clean_names(seq_len(nrow(ancillary_data)))
    )
    ancillary_data$feature_id = metabolite_data$feature_id
  } else {
    metabolite_data = ancillary_data[, character_cols]
    metabolite_data$feature_id = janitor::make_clean_names(metabolite_data[[1]])
    ancillary_data$feature_id = metabolite_data$feature_id
  }

  parsed_data$MEASUREMENTS = ancillary_data
  parsed_data$METABOLITES = metabolite_data
  return(parsed_data)
}


run_mwtab_checks = function(
  processed_mwtab,
  min_n = 5,
  min_metabolites = 100,
  min_ssf = 2
) {
  # processed_mwtab = tar_read(processed_AN002428)
  # processed_mwtab = tar_read(processed_AN003177)
  # processed_mwtab = tar_read(processed_AN003224)
  # processed_mwtab = tar_read(processed_AN003894)
  # min_n = 5
  # min_metabolites = 100
  # min_ssf = 2

  check_res = list(
    FEATURE_CHECK = "",
    SSF_CHECK = "",
    NA_CHECK = "",
    RANK_CHECK = "",
    ID = processed_mwtab$ID
  )
  if (is.null(check_res$ID)) {
    check_res$ID = ""
  }
  n_features = nrow(processed_mwtab$MEASUREMENTS)
  if (is.na(n_features) || is.null(n_features)) {
    check_res$FEATURE_CHECK = "NA FEATURES"
    processed_mwtab$CHECK = check_res
    return(processed_mwtab)
  }

  if (n_features < min_metabolites) {
    check_res$FEATURE_CHECK = "FEW METABOLITES"
    processed_mwtab$CHECK = check_res
    return(processed_mwtab)
  } else {
    check_res$FEATURE_CHECK = "GOOD"
  }

  if (is.null(processed_mwtab$METABOLITES)) {
    feature_metadata = data.frame(
      metabolite_name = processed_mwtab$MEASUREMENTS$id,
      feature_id = janitor::make_clean_names(processed_mwtab$MEASUREMENTS$id)
    )
    processed_mwtab$METABOLITES = feature_metadata
  } else {
    processed_mwtab$METABOLITES$feature_id = janitor::make_clean_names(
      processed_mwtab$METABOLITES$metabolite_name
    )
  }

  ssf_data = processed_mwtab$SUBJECT_SAMPLE_FACTORS
  match_samples = base::intersect(
    ssf_data$sample_id,
    colnames(processed_mwtab$MEASUREMENTS)
  )
  ssf_data2 = ssf_data |>
    dplyr::filter(sample_id %in% match_samples)

  ssf_check_result = check_ssf(ssf_data2, min_ssf, min_n)
  processed_mwtab$SUBJECT_SAMPLE_FACTORS = ssf_data2

  processed_mwtab$MEASUREMENTS = processed_mwtab$MEASUREMENTS[, c(
    "feature_id",
    "id",
    ssf_data2$sample_id
  )]

  check_res$SSF_CHECK = ssf_check_result
  if (ssf_check_result %in% "NOT ENOUGH SAMPLES") {
    processed_mwtab$CHECK = check_res
    return(processed_mwtab)
  }

  measurements = processed_mwtab$MEASUREMENTS[, ssf_data2$sample_id] |>
    as.matrix()
  factors = ssf_data2$factors

  missingness_rank_check = check_missingness_ranks(measurements, factors)

  check_res$RANK_CHECK = missingness_rank_check

  processed_mwtab$CHECK = check_res
  processed_mwtab
}

check_ssf = function(ssf_data2, min_ssf = 2, min_n = 5) {
  ssf_n = ssf_data2 |>
    dplyr::summarise(n = dplyr::n(), .by = factors)

  if (any(ssf_n$n >= min_n) && (nrow(ssf_n) >= min_ssf)) {
    return("GOOD SSF")
  } else {
    return("NOT ENOUGH SAMPLES")
  }
}


count_factors_replicates_mwtab = function(subject_sample_factors) {
  # subject_sample_factors = parsed_data$SUBJECT_SAMPLE_FACTORS

  other_cols = names(subject_sample_factors)[
    !(names(subject_sample_factors) %in% "sample_id")
  ]

  concat_factors = subject_sample_factors |>
    tidyr::unite("factors", tidyselect::all_of(other_cols), sep = "|")
  subject_sample_factors = dplyr::left_join(
    subject_sample_factors,
    concat_factors,
    by = "sample_id"
  )

  subject_sample_factors$sample_id = janitor::make_clean_names(
    subject_sample_factors$sample_id
  )

  n_samples = subject_sample_factors |>
    dplyr::summarise(n = dplyr::n(), .by = factors)
  list(sample_factors = subject_sample_factors, n = n_samples)
}

parse_tabs = function(in_block) {
  # in_block = block_data$STUDY
  tab_split = stringr::str_split_fixed(in_block, "\t", 2) |>
    as.data.frame() |>
    tibble::as_tibble()
  names(tab_split) = c("field", "value")

  tab_split
}

parse_nmr_data_mwtab = function(nmr_data) {
  # nmr_data = block_data[["NMR_BINNED_DATA"]]
  start_loc = which(stringr::str_detect(nmr_data, "NMR_BINNED_DATA_START")) +
    1
  end_loc = length(nmr_data) - 1
  nmr_data = nmr_data[seq(start_loc, end_loc)]

  if (length(nmr_data) <= 5) {
    return(NULL)
  }

  nmr_split = stringr::str_split(nmr_data, "\t")

  nmr_nobj = purrr::map_int(nmr_split, length)

  if (!all(nmr_nobj[1] == nmr_nobj)) {
    return(NULL)
  }
  nmr_df = purrr::map(nmr_split, \(in_split) {
    tibble::tibble(id = in_split[1], value = in_split[seq(2, nmr_nobj[1])])
  })

  is_numeric = purrr::map_lgl(nmr_df, \(in_df) {
    num_vals = suppressWarnings(as.numeric(in_df[["value"]]))
    sum(!is.na(num_vals)) > 1
  })

  samples_df = which(purrr::map_lgl(nmr_df[1:5], \(in_df) {
    grepl("Bin range|Samples", in_df[["id"]][1], ignore.case = TRUE)
  }))

  is_numeric[samples_df] = FALSE

  nmr_numeric = nmr_df[is_numeric]

  nmr_long_df = purrr::map(nmr_numeric, \(in_numeric) {
    in_numeric$sample_id = nmr_df[[samples_df]]$value
    in_numeric
  }) |>
    purrr::list_rbind()
  nmr_long_df$value = as.numeric(nmr_long_df$value)
  nmr_wide_df = nmr_long_df |>
    tidyr::pivot_wider(
      id_cols = id,
      values_from = value,
      names_from = sample_id
    )
  names(nmr_wide_df) = janitor::make_clean_names(names(nmr_wide_df))
  nmr_wide_df$feature_id = janitor::make_clean_names(nmr_wide_df$id)
  nmr_wide_df
}

parse_ms_data_mwtab = function(ms_data) {
  # ms_data = block_data[["MS_METABOLITE_DATA"]]
  # message("ms data!")
  start_loc = which(stringr::str_detect(ms_data, "MS_METABOLITE_DATA_START")) +
    1
  end_loc = length(ms_data) - 1

  ms_data = ms_data[seq(start_loc, end_loc)]

  if (length(ms_data) <= 5) {
    return(NULL)
  }
  ms_split = stringr::str_split(ms_data, "\t")

  ms_nobj = purrr::map_int(ms_split, length)

  if (!all(ms_nobj[1] == ms_nobj)) {
    return(NULL)
  }
  ms_df = purrr::map(ms_split, \(in_split) {
    tibble::tibble(id = in_split[1], value = in_split[seq(2, ms_nobj[1])])
  })

  is_numeric = purrr::map_lgl(ms_df, \(in_df) {
    num_vals = suppressWarnings(as.numeric(in_df[["value"]]))
    sum(!is.na(num_vals)) > 0
  })

  if (sum(is_numeric) < 2) {
    return(NULL)
  }

  samples_df = which(purrr::map_lgl(ms_df[1:5], \(in_df) {
    in_df[["id"]][1] == "Samples"
  }))

  if (any(nchar(ms_df[[samples_df]]$value) == 0)) {
    return(NULL)
  }
  is_numeric[samples_df] = FALSE

  ms_numeric = ms_df[is_numeric]

  ms_long_df = purrr::map(ms_numeric, \(in_numeric) {
    in_numeric$sample_id = ms_df[[samples_df]]$value
    in_numeric
  }) |>
    purrr::list_rbind()
  ms_long_df$value = as.numeric(ms_long_df$value)

  # sometimes there are multiple entries
  mult_entry = ms_long_df |>
    dplyr::summarise(n = dplyr::n(), .by = c(id, sample_id)) |>
    dplyr::filter(n > 1L) |>
    dplyr::pull(id) |>
    unique()
  ms_long_df_uniq = ms_long_df |>
    dplyr::filter(!(id %in% mult_entry))

  if (nrow(ms_long_df_uniq) == 0) {
    return(NULL)
  }

  ms_wide_df = ms_long_df_uniq |>
    tidyr::pivot_wider(
      id_cols = id,
      values_from = value,
      names_from = sample_id
    )

  names(ms_wide_df) = janitor::make_clean_names(names(ms_wide_df))
  ms_wide_df$feature_id = janitor::make_clean_names(ms_wide_df$id)

  ms_wide_df
}


parse_factors_mwtab = function(factor_data) {
  # message("factor data!")
  # factor_data = block_data[["SUBJECT_SAMPLE_FACTORS"]]

  factor_out = stringr::str_split(factor_data, "\t")
  factor_data = purrr::map(factor_out, \(in_factor) {
    sub_sample = tibble::tibble(
      subject = in_factor[2],
      sample_id = in_factor[3]
    )
    factors = split_and_return_factors(in_factor[4])
    out_df = dplyr::bind_cols(sub_sample, factors)
    out_df
  }) |>
    purrr::list_rbind()
  null_subject = factor_data[["subject"]] == "-"
  if (any(factor_data$factor %in% "sample_id")) {
    return(NULL)
  }

  # these basically can't be empty, if they are, bad things happen.
  if (any(nchar(factor_data$factor) == 0)) {
    return(NULL)
  }
  if (any(nchar(factor_data$sample_id) == 0)) {
    return(NULL)
  }
  if (sum(null_subject) == nrow(factor_data)) {
    factor_data$subject = NULL
    wider_factors = factor_data |>
      tidyr::pivot_wider(
        id_cols = sample_id,
        names_from = factor,
        values_from = item
      )
  } else {
    wider_factors = factor_data |>
      tidyr::pivot_wider(
        id_cols = c(subject, sample_id),
        names_from = factor,
        values_from = item
      )
  }

  return(wider_factors)
}

split_and_return_factors = function(factor_string) {
  # factor_string = "AFTER_MEAL_TIME:30_MIN | gender:female | POST_SURGERY_TIME:12MO | SURGERY:GASTRIC_BAND"
  # factor_string = "Treatment:Media alone"

  pipe_split = stringr::str_split(factor_string, "\\|")[[1]] |> trimws()
  colon_split = stringr::str_split(pipe_split, "\\:")
  factor_df = tibble::tibble(
    factor = purrr::map_chr(colon_split, \(x) {
      x[1]
    }),
    item = purrr::map_chr(colon_split, \(x) {
      x[2]
    })
  )
  factor_df
}


parse_metabolites_mwtab = function(metabolites_block) {
  # metabolites_block = block_data$METABOLITES

  end_loc = which(grepl("METABOLITES_END", metabolites_block)) - 1
  metabolites_range = seq(2, end_loc)
  metabolites_block = metabolites_block[metabolites_range]

  if (length(metabolites_block) < 2) {
    return(NULL)
  }
  has_id_row = grepl(
    "metabolite_name|metabolite_id|id|kegg",
    metabolites_block[1],
    ignore.case = TRUE
  )

  all_split = stringr::str_split(metabolites_block, "\t")
  if (has_id_row) {
    col_id = janitor::make_clean_names(all_split[[1]])
    all_split = all_split[-1]
  } else {
    col_id = janitor::make_clean_names(seq(1, length(all_split[[1]])))
  }

  split_bind = do.call(rbind, all_split)

  if (ncol(split_bind) < length(col_id)) {
    col_id = col_id[seq_len(ncol(split_bind))]
  } else if (ncol(split_bind) > length(col_id)) {
    col_id = c(
      col_id,
      janitor::make_clean_names(seq_len(ncol(split_bind) - length(col_id)))
    )
  }
  colnames(split_bind) = col_id
  tibble::as_tibble(split_bind)
}


convert_json_smd = function(rds_file, reps, smd_file, ...) {
  # rep_data = readRDS("data/processed/mwtab_keep.rds")
  # use_reps = rep_data |>
  #   dplyr::filter(id %in% "AN000023") |>
  #   dplyr::pull(reps)
  # use_reps = use_reps[[1]]
  # mwtab_data = readRDS("data/processed/AN000023.rds")
  # mwtab_keep = readRDS("data/processed/mwtab_keep.rds")
  # rds_file = mwtab_keep$rds_file[24]
  # reps = mwtab_keep$reps[[24]]
  # smd_file = mwtab_keep$smd_file[24]

  # rds_file = mwtab_keep$rds_file[6]
  # reps = mwtab_keep$reps[[6]]
  # smd_file = mwtab_keep$smd_file[6]

  # rds_file = mwtab_keep$rds_file[7]
  # reps = mwtab_keep$reps[[7]]
  # smd_file = mwtab_keep$smd_file[7]

  # rds_file = mwtab_keep$rds_file[118]
  # reps = mwtab_keep$reps[[118]]
  # smd_file = mwtab_keep$smd_file[118]
  mwtab_data = readRDS(rds_file)
  sample_metadata = mwtab_data$SUBJECT_SAMPLE_FACTORS
  sample_metadata = sample_metadata |>
    dplyr::filter(factors %in% reps$factors)
  if (mwtab_data$MEASUREMENT_TYPE %in% "MS") {
    mwtab_data$MEASUREMENT_INFO = mwtab_data$MS
  } else if (mwtab_data$MEASUREMENT_TYPE %in% "NMR") {
    mwtab_data$MEASUREMENT_INFO = mwtab_data$NMR
  }
  project_data = mwtab_data[c(
    "MWINFO",
    "PROJECT",
    "STUDY",
    "SUBJECT",
    "COLLECTION",
    "TREATMENT",
    "SAMPLEPREP",
    "CHROMATOGRAPHY",
    "ANALYSIS",
    "METABOLITES",
    "N_REPLICATES",
    "MEASUREMENTS",
    "MEASUREMENT_TYPE",
    "MEASUREMENT_INFO"
  )]

  if (!is.null(mwtab_data$METABOLITES)) {
    feature_metadata = mwtab_data$METABOLITES
  } else {
    feature_metadata = data.frame(
      metabolite_name = mwtab_data$MEASUREMENTS$id
    )
  }

  match_samples = base::intersect(
    sample_metadata$sample_id,
    colnames(mwtab_data$MEASUREMENTS)
  )

  if ((length(match_samples) == 0) || (is.null(match_samples))) {
    stop("No matching samples!")
  }
  sample_metadata = sample_metadata |>
    dplyr::filter(sample_id %in% match_samples)
  measurements = mwtab_data$MEASUREMENTS[, c("id", sample_metadata$sample_id)]

  feature_DF = DataFrame(feature_metadata)
  feature_DF$feature_id = janitor::make_clean_names(feature_DF$metabolite_name)
  rownames(feature_DF) = feature_DF$feature_id

  sample_DF = DataFrame(sample_metadata)
  sample_DF$sample_id = janitor::make_clean_names(sample_DF$sample_id)
  colnames(sample_DF) = janitor::make_clean_names(colnames(sample_DF))
  rownames(sample_DF) = sample_DF$sample_id
  measurement_array = measurements |>
    dplyr::select(-id) |>
    as.matrix()
  colnames(measurement_array) = janitor::make_clean_names(colnames(
    measurement_array
  ))
  rownames(measurement_array) = janitor::make_clean_names(measurements$id)
  keep_metabolites = base::intersect(
    rownames(feature_DF),
    rownames(measurement_array)
  )
  feature_DF = feature_DF[keep_metabolites, ]
  measurement_array = measurement_array[
    feature_DF$feature_id,
    sample_DF$sample_id
  ]

  measurement_array[measurement_array <= 0] = NA

  out_smd = SummarizedExperiment(
    assays = SimpleList(counts = measurement_array),
    rowData = feature_DF,
    colData = sample_DF,
    metadata = project_data
  )
  saveRDS(out_smd, smd_file)
}


convert_nsclc_smd = function() {
  data_in = readRDS("data/nsclc/nsclc_count_info.rds")
  lipid_in = readRDS("data/nsclc/nsclc_feature_lipid_classes.rds")

  sample_DF = data_in$info |>
    dplyr::mutate(
      sample_id = janitor::make_clean_names(sample),
      factors = treatment
    ) |>
    janitor::clean_names() |>
    DataFrame()
  rownames(sample_DF) = sample_DF$sample_id

  measurements_array = data_in$counts
  colnames(measurements_array) = janitor::make_clean_names(colnames(
    measurements_array
  ))
  rownames(measurements_array) = janitor::make_clean_names(rownames(
    measurements_array
  ))

  measurements_features = data.frame(feature_id = rownames(measurements_array))
  lipid_distinct = dplyr::distinct(lipid_in)
  feature_df = dplyr::full_join(
    lipid_distinct,
    measurements_features,
    by = "feature_id"
  )
  testthat::expect_equal(nrow(feature_df), nrow(measurements_array))
  feature_DF = feature_df |> janitor::clean_names() |> DataFrame()
  rownames(feature_DF) = feature_DF$feature_id

  measurements_array = measurements_array[
    feature_DF$feature_id,
    sample_DF$sample_id
  ]

  measurements_array[measurements_array == 0] = NA

  out_smd = SummarizedExperiment(
    assays = SimpleList(counts = measurements_array),
    rowData = feature_DF,
    colData = sample_DF
  )
  saveRDS(out_smd, "data/smd/NSCLC.rds")

  mwtab_keep = readRDS("data/smd/mwtab_smd.rds")
  mwtab_keep = dplyr::bind_rows(
    mwtab_keep,
    tibble::tibble(
      id = "NSCLC",
      type = "MS",
      reps = NULL,
      n_features = nrow(out_smd),
      rds_file = "",
      value_check = "good",
      smd_file = "data/smd/NSCLC.rds"
    )
  )
  saveRDS(mwtab_keep, "data/smd/mwtab_smd.rds")
}
