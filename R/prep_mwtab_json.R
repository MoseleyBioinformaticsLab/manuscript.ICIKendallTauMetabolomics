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


parse_json = function(mwtab_file, id, ancillary_path) {
  #mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000023%2Fmwtab%2Ftxt.txt"  # MS file
  # mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000693%2Fmwtab%2Ftxt.txt"  # NMR file

  # mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000001%2Fmwtab%2Ftxt.txt"

  #mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000004%2Fmwtab%2Ftxt.txt"

  # mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000038%2Fmwtab%2Ftxt.txt"

  # mwtab_file = "data/repaired/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000555%2Fmwtab%2Ftxt.txt"

  # mwtab_file = tar_read(dataset_AN000172)
  # mwtab_file = tar_read(mwtab_AN003177)
  # mwtab_file = tar_read(mwtab_AN002445)
  # mwtab_file = tar_read(mwtab_AN003894)
  # mwtab_file = tar_read(mwtab_AN001238)
  # mwtab_file = tar_read(dataset_AN000023)
  # mwtab_file = tar_read(dataset_AN000555)
  # mwtab_file = tar_read(dataset_AN001154)
  # mwtab_file = tar_read(dataset_AN007151)
  # id = "AN007151"
  # mwtab_file = tar_read(dataset_AN004150)
  # id = "AN004150"
  # mwtab_file = tar_read(dataset_AN003938)
  # id = "AN003938"
  # mwtab_file = tar_read(dataset_AN000428)
  # id = "AN000428"
  # mwtab_file = tar_read(dataset_AN000111)
  # id = "AN000111"
  # mwtab_file = tar_read(dataset_AN000352) # check on this one to see what is going on
  # id = "AN000352"
  # mwtab_file = tar_read(dataset_AN000593)
  # id = "AN000593"
  mwtab_list = jsonlite::fromJSON(mwtab_file, simplifyVector = FALSE)

  parsed_data = purrr::imap(mwtab_list, \(block, block_id) {
    switch(
      block_id,
      `METABOLOMICS WORKBENCH` = parse_list(block),
      PROJECT = parse_list(block),
      STUDY = parse_list(block),
      SUBJECT = parse_list(block),
      SUBJECT_SAMPLE_FACTORS = parse_factors_json(block),
      COLLECTION = parse_list(block),
      TREATMENT = parse_list(block),
      SAMPLEPREP = parse_list(block),
      CHROMATOGRAPHY = parse_list(block),
      ANALYSIS = parse_list(block),
      MS = parse_list(block),
      MS_METABOLITE_DATA = suppressMessages(parse_ms_data_json(block)),
      NM = parse_list(block),
      NMR_BINNED_DATA = suppressMessages(parse_nmr_data_json(block)),
      NMR_METABOLITE_DATA = suppressMessages(parse_nmr_data_json(block)),
      block
    )
  })

  parsed_data$ID = id
  if (is.null(parsed_data$SUBJECT_SAMPLE_FACTORS)) {
    # browser()
    return(parsed_data)
  }

  if (!is.null(parsed_data$MS_METABOLITE_DATA)) {
    parsed_data$MEASUREMENTS = parsed_data$MS_METABOLITE_DATA$DATA
    parsed_data$METABOLITES = parsed_data$MS_METABOLITE_DATA$METABOLITES
    parsed_data$MEASUREMENT_TYPE = "MS"
  } else if (!is.null(parsed_data$NMR_BINNED_DATA)) {
    parsed_data$MEASUREMENTS = parsed_data$NMR_BINNED_DATA$DATA
    parsed_data$METABOLITES = parsed_data$NMR_BINNED_DATA$METABOLITES
    parsed_data$MEASUREMENT_TYPE = "NMR"
  } else if (!is.null(parsed_data$NMR_METABOLITE_DATA)) {
    parsed_data$MEASUREMENTS = parsed_data$NMR_METABOLITE_DATA$DATA
    parsed_data$METABOLITES = parsed_data$NMR_METABOLITE_DATA$METABOLITES
    parsed_data$MEASUREMENT_TYPE = "NMR"
  } else {
    parsed_data = check_and_parse_ancillary_json(parsed_data, ancillary_path)
  }

  parsed_data
}

parse_list = function(in_list) {
  tibble::tibble(field = names(in_list), value = unlist(in_list))
}

create_results_txt = function(mwtab_df) {
  study_id = mwtab_df |>
    dplyr::filter(field %in% "STUDY_ID") |>
    dplyr::pull(value)

  analysis_id = mwtab_df |>
    dplyr::filter(field %in% "ANALYSIS_ID") |>
    dplyr::pull(value)

  if ((length(study_id) == 0) || (length(analysis_id) == 0)) {
    return(character(0))
  }
  results_txt = paste0(study_id, "_", analysis_id, "_Results.txt")
  return(results_txt)
}

check_and_parse_ancillary_json = function(
  parsed_data,
  ancillary_path = "/big_data/data/MetabolomicsWorkbench/ancillary/"
) {
  # parsed_data = tar_read(processed_AN000111)
  # ancillary_path = "/big_data/data/rmflight_icikt_poster_iecm/data/ancillary/"
  if (!is.null(parsed_data$MS)) {
    parsed_data$MEASUREMENT_TYPE = "MS"
  }
  if (!is.null(parsed_data$NM)) {
    parsed_data$MEASUREMENT_TYPE = "NMR"
  }
  results_file = character(0)

  results_file = create_results_txt(parsed_data[["METABOLOMICS WORKBENCH"]])

  if (length(results_file) == 0) {
    return(parsed_data)
  }

  ancillary_files = fs::dir_ls(ancillary_path)
  ancillary_just_file = fs::path_file(ancillary_files)
  names(ancillary_files) = ancillary_just_file

  match_file = ancillary_files[results_file]

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

  check_character = function(data_df) {
    col_types = readr::spec(data_df)
    is_character = purrr::map_lgl(col_types$cols, \(x) {
      any(grepl("character", class(x)))
    })
    is_character
  }
  check_character_retype = function(ancillary_data) {
    character_cols = check_character(ancillary_data)

    if (sum(character_cols) == ncol(ancillary_data)) {
      ancillary_tmp = ancillary_data[-1, ] |>
        suppressMessages(readr::type_convert())
      character_cols = check_character(ancillary_tmp)
      if (sum(character_cols) == ncol(ancillary_tmp)) {
        ancillary_tmp = ancillary_data[c(-1, -2), ] |>
          suppressMessages(readr::type_convert())
        character_cols = check_character(ancillary_tmp)
        if (sum(character_cols) == ncol(ancillary_data)) {
          return(NULL)
        } else {
          return(ancillary_tmp)
        }
      } else {
        return(ancillary_tmp)
      }
    } else {
      return(ancillary_data)
    }
  }
  # not everything can be a character column, there should
  # be some numeric columns

  # we assume the first character column is an ID.
  # if there are no character columns, we just make our own.
  ancillary_checked = check_character_retype(ancillary_data)
  if (is.null(ancillary_checked)) {
    return(parsed_data)
  }
  is_character = check_character(ancillary_checked)
  if (sum(is_character) == 0) {
    metabolite_data = tibble::tibble(
      feature_id = janitor::make_clean_names(seq_len(nrow(ancillary_checked)))
    )
    ancillary_checked$feature_id = metabolite_data$feature_id
  } else {
    metabolite_data = ancillary_checked[, is_character]
    metabolite_data$feature_id = janitor::make_clean_names(metabolite_data[[1]])
    ancillary_checked$feature_id = metabolite_data$feature_id
  }

  parsed_data$MEASUREMENTS = ancillary_checked
  parsed_data$METABOLITES = metabolite_data
  return(parsed_data)
}


get_check = function(checked_mwtab) {
  tibble::as_tibble(checked_mwtab$CHECK)
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

check_missingness_ranks = function(measurements, factors) {
  if ((sum(is.na(measurements)) == 0) && (sum(measurements == 0) == 0)) {
    return("NO MISSING VALUES")
  }

  ranked_data = ICIKendallTau::rank_order_data(
    measurements,
    sample_classes = factors
  )

  ranked_null = purrr::map_lgl(ranked_data, is.null)
  if (any(ranked_null)) {
    return("ALL NA IN A GROUP")
  }

  ranked_cor = purrr::map(ranked_data, \(in_data) {
    min_rank = in_data$n_na_rank |>
      dplyr::summarise(min_rank = min(median_rank), .by = n_na)
    full_cor = cor(
      in_data$n_na_rank$n_na,
      in_data$n_na_rank$median_rank,
      method = "kendall"
    )
    min_cor = cor(
      min_rank$n_na,
      min_rank$min_rank,
      method = "kendall"
    )
    tibble::tibble(median = full_cor, min = min_cor)
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(sign_diff = sign(median) * sign(min))

  if (any(is.na(ranked_cor$sign_diff))) {
    return("NA CORRELATION")
  } else if (any(ranked_cor$sign_diff == -1)) {
    return("SIGN DIFFERENCE")
  } else if (any(sign(ranked_cor$median) == 1)) {
    return("POSITIVE MEDIAN")
  } else if (any(sign(ranked_cor$min) == 1)) {
    return("POSITIVE MINIMUM")
  } else {
    return("GOOD")
  }
}

count_factors_replicates_long = function(subject_sample_factors) {
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


parse_nmr_data_json = function(nmr_data) {
  # nmr_data = block_data[["NMR_BINNED_DATA"]]
  # nmr_data = mwtab_list$NMR_BINNED_DATA

  nmr_df_1 = try(tibble::as_tibble(nmr_data$Data[[1]]))
  if (inherits(nmr_df_1, "try-error")) {
    return(NULL)
  }

  nmr_data_df = tryCatch(
    {
      purrr::map(nmr_data$Data, tibble::as_tibble) |>
        purrr::list_rbind() |>
        readr::type_convert() |>
        janitor::clean_names()
    },
    error = return(NULL)
  )

  nmr_data_df$feature_id = janitor::make_clean_names(nmr_data_df[[1]])

  if (!is.null(nmr_data$Metabolites)) {
    metabolite_df = tryCatch(
      {
        purrr::map(nmr_data$Metabolites, tibble::as_tibble) |>
          purrr::list_rbind() |>
          readr::type_convert() |>
          janitor::clean_names()
      },
      error = return(NULL)
    )
  } else {
    metabolite_df = tibble::tibble(feature_id = nmr_data_df$feature_id)
  }

  list(DATA = nmr_data_df, METABOLITES = metabolite_df)
}

parse_ms_data_json = function(ms_data) {
  # ms_data = block_data[["MS_METABOLITE_DATA"]]
  # message("ms data!")
  # ms_data = mwtab_list$MS_METABOLITE_DATA

  ms_df_1 = try(tibble::as_tibble(ms_data$Data[[1]]))
  if (inherits(ms_df_1, "try-error")) {
    return(NULL)
  }
  ms_data_df = tryCatch(
    {
      purrr::map(ms_data$Data, tibble::as_tibble) |>
        purrr::list_rbind() |>
        readr::type_convert() |>
        janitor::clean_names()
    },
    error = return(NULL)
  )
  ms_data_df$feature_id = janitor::make_clean_names(ms_data_df[[1]])

  if (!is.null(ms_data$Metabolites)) {
    metabolite_df = tryCatch(
      {
        purrr::map(ms_data$Metabolites, tibble::as_tibble) |>
          purrr::list_rbind() |>
          readr::type_convert() |>
          janitor::clean_names()
      },
      error = return(NULL)
    )
  } else {
    metabolite_df = tibble::tibble(feature_id = ms_data_df$feature_id)
  }

  return(list(DATA = ms_data_df, METABOLITES = metabolite_df))
}


parse_factors_json = function(factor_data) {
  # message("factor data!")
  # factor_data = mwtab_list[["SUBJECT_SAMPLE_FACTORS"]]

  factor_df = purrr::map(factor_data, \(in_sample) {
    # in_sample = factor_data[[1]]
    sub_sample = tibble::tibble(
      subject = in_sample[["Subject ID"]],
      sample_id = janitor::make_clean_names(in_sample[["Sample ID"]])
    )
    factors = tibble::tibble(
      factor = names(in_sample$Factors),
      item = unlist(in_sample$Factors),
      type = "subject_sample_factor"
    )

    out_df = dplyr::bind_cols(sub_sample, factors)

    if (!is.null(in_sample$`Additional sample data`)) {
      additional = tibble::tibble(
        factor = names(in_sample$`Additional sample data`),
        item = unlist(in_sample$`Additional sample data`),
        type = "additional_sample_data"
      )
      add_df = dplyr::bind_cols(sub_sample, additional)
      out_df = dplyr::bind_rows(out_df, add_df)
    }
    out_df
  }) |>
    purrr::list_rbind()

  null_subject = factor_df[["subject"]] == "-"
  if (any(factor_df$factor %in% "sample_id")) {
    return(NULL)
  }

  # these basically can't be empty, if they are, bad things happen.
  if (any(nchar(factor_df$factor) == 0)) {
    return(NULL)
  }
  if (any(nchar(factor_df$sample_id) == 0)) {
    return(NULL)
  }

  return(factor_df)
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


parse_metabolites_block = function(metabolites_block) {
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


find_possible_mwtab = function() {
  mwtab_results = readRDS("data/processed/mwtab_results.rds")
  n_analyses = nrow(mwtab_results)
  n_na = sum(is.na(mwtab_results$n_features))
  n_less_200 = sum(
    !is.na(mwtab_results$n_features) & (mwtab_results$n_features < 200),
    na.rm = TRUE
  )
  mwtab_nonna = mwtab_results |>
    dplyr::filter(!is.na(n_features)) |>
    dplyr::filter(n_features >= 200)
  n_nonna = nrow(mwtab_nonna)
  keep_reps = purrr::map(mwtab_nonna$reps, \(in_reps) {
    # in_reps = mwtab_nonna$reps[[1]]
    keep_reps = in_reps |>
      dplyr::filter(n >= 3)
    keep_reps
  })
  has_multi = purrr::map_lgl(keep_reps, \(x) {
    (nrow(x) >= 2) && (sum(x$n >= 5) >= 1)
  })
  n_has_multi = sum(has_multi)

  mwtab_nonna$reps = keep_reps

  check_results = mwtab_nonna[has_multi, ]
  n_data = tibble::tibble(
    `N Analyses` = n_analyses,
    `N NA Features` = n_na,
    `N Less 200` = n_less_200,
    `N Few Factors` = sum(!has_multi),
    `N Multiple Factors` = n_has_multi
  )
  out_results = list(N = n_data, check = check_results)
  saveRDS(out_results, "data/processed/mwtab_check.rds")
}


check_ranking_mwtab = function() {
  out_results = readRDS("data/processed/mwtab_check.rds")

  check_results = out_results$check |>
    dplyr::mutate(rds_file = fs::path("data", "processed", id, ext = "rds"))

  check_each = purrr::pmap_chr(
    check_results,
    load_and_check_mwtab,
    .progress = TRUE
  )

  check_results$value_check = check_each
  saveRDS(check_results, "data/processed/mwtab_rank.rds")
  return(invisible(NULL))
}

find_mwtab_keep = function() {
  check_results = readRDS("data/processed/mwtab_rank.rds")

  keep_results = check_results |>
    dplyr::filter(value_check %in% c("good", "sign difference"))

  saveRDS(keep_results, "data/processed/mwtab_keep.rds")
}

load_and_check_mwtab = function(rds_file, reps, ...) {
  # rds_file = check_results$rds_file[70]
  # reps = check_results$reps[[70]]

  # rds_file = check_results$rds_file[7]
  # reps = check_results$reps[[7]]

  keep_factors = reps
  mwtab_data = readRDS(rds_file)

  sample_factors = mwtab_data$SUBJECT_SAMPLE_FACTORS |>
    dplyr::filter(factors %in% keep_factors$factors)

  sample_subjects = mwtab_data$SUBJECT
  if (any(grepl("\\;", sample_subjects$value))) {
    return("multiple species")
  }

  if (
    length(
      base::setdiff(
        sample_factors$sample_id,
        colnames(mwtab_data$MEASUREMENTS)
      ) >
        0
    )
  ) {
    match_samples = base::intersect(
      sample_factors$sample_id,
      colnames(mwtab_data$MEAUREMENTS)
    )
    sample_factors = sample_factors |>
      dplyr::filter(sample_id %in% match_samples)

    n_rep = sample_factors |>
      dplyr::summarise(n = dplyr::n(), .by = factors)
    keep_it = (nrow(n_rep) >= 2) && (sum(n_rep$n >= 5) >= 1)

    if (!keep_it) {
      return("missing or not enough samples")
    }
  }

  if (
    inherits(mwtab_data$MEASUREMENTS[[sample_factors$sample_id[1]]], "list")
  ) {
    return("list data")
  }

  metabolite_values = mwtab_data$MEASUREMENTS[, sample_factors$sample_id] |>
    as.matrix()

  ranked_data = ICIKendallTau::rank_order_data(
    metabolite_values,
    sample_classes = sample_factors$factors
  )
  ranked_null = purrr::map_lgl(ranked_data, is.null)
  if (any(ranked_null)) {
    return("all NA in a group")
  }

  if (
    (sum(is.na(metabolite_values)) == 0) && (sum(metabolite_values == 0) == 0)
  ) {
    return("no missing")
  }
  rank_cor = purrr::map(ranked_data, \(in_data) {
    # in_data = ranked_data[[1]]
    min_rank = in_data$n_na_rank |>
      dplyr::summarise(min_rank = min(median_rank), .by = n_na)
    full_cor = cor(
      in_data$n_na_rank$n_na,
      in_data$n_na_rank$median_rank,
      method = "kendall"
    )
    min_cor = cor(
      min_rank$n_na,
      min_rank$min_rank,
      method = "kendall"
    )
    tibble::tibble(median = full_cor, min = min_cor)
  }) |>
    purrr::list_rbind()

  rank_cor = rank_cor |>
    dplyr::mutate(sign_diff = sign(median) * sign(min))

  if (any(is.na(rank_cor$sign_diff))) {
    return("na corr")
  } else if (any(rank_cor$sign_diff == -1)) {
    return("sign difference")
  } else if (any(sign(rank_cor$median) == 1)) {
    return("positive median")
  } else if (any(sign(rank_cor$min) == 1)) {
    return("positive minimum")
  } else {
    return("good")
  }
}


convert_mwtab_json_smd = function(rds_file, reps, smd_file, ...) {
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


convert_list_to_df = function(measurements) {
  id = "AN000555"
  rds_file = fs::path("data", "processed", id, ext = "rds")
  mwtab_data = readRDS(rds_file)
  measurements = mwtab_data$MEASUREMENTS

  out_vals = purrr::map(measurements, \(in_meas) {
    unlist(in_meas)
  })
}


create_ancillary_file_list = function(repaired_json_dir) {
  repaired_json_dir = "/big_data/data/MetabolomicsWorkbench/repaired/json_2025-11-12/"
  json_files = fs::dir_ls(repaired_json_dir)
  tar_source("R")

  check_external_results = function(in_file) {
    # in_file = "/big_data/data/MetabolomicsWorkbench/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN003030%2Fmwtab%2Ftxt.txt.json"
    # in_file = "/big_data/data/MetabolomicsWorkbench/repaired/json_2025-11-12/https:%2F%2Fwww_metabolomicsworkbench_org%2Frest%2Fstudy%2Fanalysis_id%2FAN000001%2Fmwtab%2Ftxt.txt.json"
    json_data = jsonlite::fromJSON(in_file, simplifyVector = FALSE)
    results_file = paste0(
      json_data$`METABOLOMICS WORKBENCH`$STUDY_ID,
      "_",
      json_data$`METABOLOMICS WORKBENCH`$ANALYSIS_ID,
      "_Results.txt"
    )
    tmp_id = extract_mwtab_ids(in_file)
    tmp_id$results = results_file
    tmp_id
  }

  results_files = purrr::map(
    json_files,
    check_external_results,
    .progress = TRUE
  ) |>
    purrr::list_rbind()
  cat(
    results_files,
    sep = "\n",
    file = "/big_data/data/MetabolomicsWorkbench/ancillary/ancillary_files.txt"
  )
  return(invisible(NULL))
}

download_ancillary = function(
  ancillary_files = "/big_data/data/MetabolomicsWorkbench/ancillary/ancillary_files.txt",
  download_dir = "/big_data/data/MetabolomicsWorkbench/ancillary/"
) {
  # the list of files can be created by doing:
  # grep -R -h -o "ST.*[[:digit:]].*\\_AN.*[[:digit:]].*.txt" > ancillary_files.txt
  url_loc = "https://www.metabolomicsworkbench.org/studydownload/"
  exist_files = fs::dir_ls(download_dir) |> fs::path_file()
  all_files = readLines(ancillary_files)
  all_files = all_files[!(all_files %in% exist_files)]
  url_full = paste0(url_loc, all_files)

  dest_files = fs::path(download_dir, all_files)
  for (iloc in seq_len(length(all_files))) {
    try(download.file(url = url_full[iloc], destfile = dest_files[iloc]))
    Sys.sleep(2)
  }
  return(invisible(NULL))
}


test_json_parsing = function() {
  tar_source("R")
  json_files = fs::dir_ls(
    "/big_data/data/MetabolomicsWorkbench/repaired/json_2025-11-12"
  )
  ancillary_path = "/big_data/data/MetabolomicsWorkbench/ancillary/"
  mwtab_datasets = extract_mwtab_ids(json_files)

  set.seed(1234)
  test_datasets = mwtab_datasets |>
    dplyr::slice_sample(prop = 0.025) |>
    dplyr::transmute(mwtab_file = file, id = id) |>
    dplyr::mutate(ancillary_path = ancillary_path)

  json_parsed = purrr::pmap(test_datasets, parse_json, .progress = TRUE) |>
    flighttools::ft_notify_success_error()
}
