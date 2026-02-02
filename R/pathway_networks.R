calculate_qratio = function(network, annotations, use_weights = TRUE) {
  # network = network_correlations
  # annotations = use_annotation
  # assume the network is a data.frame edges
  # sum over annotations
  # sum the entire network of weights

  # assume only the positive partial correlations are useful
  network = network |>
    dplyr::filter(weight > 0)

  if (!use_weights) {
    network = network |>
      dplyr::mutate(weight = 1)
  }
  all_features_annotations = unique(unlist(annotations))

  network_all_nodes = unique(unlist(network[, c("start_node", "end_node")]))

  network_annotated = network |>
    dplyr::filter(
      start_node %in% all_features_annotations,
      end_node %in% all_features_annotations
    )

  # network_annotated = network_annotated |>
  #   dplyr::mutate(weight = 1)

  annotated_sum = sum(network_annotated$weight)
  partitions = purrr::imap(annotations, \(in_features, id) {
    other_features = all_features_annotations[
      !(all_features_annotations %in% in_features)
    ]
    # generate self-self network - sum correlations
    within_sum = network_annotated |>
      dplyr::filter(
        ((start_node %in% in_features) & (end_node %in% in_features)) |
          ((end_node %in% in_features) & (start_node %in% in_features))
      ) |>
      dplyr::summarise(sum = sum(weight)) |>
      dplyr::pull(sum)
    out_sum = network_annotated |>
      dplyr::filter(
        ((start_node %in% in_features) & (end_node %in% other_features)) |
          (end_node %in% in_features) & (start_node %in% other_features)
      ) |>
      dplyr::summarise(sum = sum(weight)) |>
      dplyr::pull(sum)

    out_val = (within_sum / annotated_sum) - ((out_sum / annotated_sum)^2)

    n_annot = length(in_features)

    tibble::tibble(id = id, ratio = out_val, n_features = length(in_features))
  }) |>
    dplyr::bind_rows()
  q_value = sum(partitions$ratio)
  list(partitions = partitions, q_value = q_value)
}

get_just_qratio = function(feature_qratio) {
  # feature_qratio = tar_read(feature_qratio_AN005135)
  if (is.null(feature_qratio)) {
    return(NULL)
  }

  qratio_data = feature_qratio$qvalues

  out_data = purrr::imap(qratio_data, \(data, id) {
    tibble::tibble(q_value = data$q_value, correlation = id)
  }) |>
    purrr::list_rbind()
  out_data$id = feature_qratio$metadata$CHECK$ID
  out_data
}

calculate_all_network_qratios = function(
  feature_partial_correlation,
  grouped_annotations
) {
  # feature_partial_correlation = tar_read(feature_partial_correlation_AN002783)
  # grouped_annotations = tar_read(predicted_annotations_grouped)
  if (is.null(feature_partial_correlation)) {
    return(NULL)
  }

  feature_annotations = grouped_annotations |>
    dplyr::filter(id %in% feature_partial_correlation$metadata$CHECK$ID)
  tmp_features = feature_annotations |>
    dplyr::select(features) |>
    dplyr::distinct() |>
    dplyr::mutate(feature2 = janitor::make_clean_names(features))

  feature_annotations = dplyr::left_join(
    feature_annotations,
    tmp_features,
    by = "features"
  )
  use_annotations = split(
    feature_annotations$feature2,
    feature_annotations$group
  )

  prep_partial_cor = function(partial_cor_df) {
    # partial_cor_df = feature_partial_correlation$partial_cor$icikt
    partial_cor_run = partial_cor_df |>
      dplyr::filter(significant) |>
      dplyr::transmute(start_node = s1, end_node = s2, weight = partial_cor)
    partial_cor_run
  }
  qratios_all = purrr::map(feature_partial_correlation$partial_cor, \(x) {
    # x = feature_partial_correlation$partial_cor[[1]]
    use_cor = prep_partial_cor(x)
    calculate_qratio(use_cor, use_annotations)
  })
  feature_partial_correlation$qvalues = qratios_all
  feature_partial_correlation
}

calculate_feature_network_qratio_new = function(
  partial_correlations,
  annotations,
  use_weights = TRUE
) {
  if (is.null(partial_correlations)) {
    return(NULL)
  }
  network_correlations = partial_correlations |>
    dplyr::filter(significant) |>
    dplyr::transmute(start_node = s1, end_node = s2, weight = partial_cor)
}

calculate_feature_network_qratio = function(
  partial_correlations,
  annotations,
  compound_type = "pathway",
  compound_mapping = metabolite_kegg,
  lipid_type = "categories",
  use_weights = TRUE
) {
  # partial_correlations = tar_read(feature_partial_cor_ici_yeast)
  # partial_correlations = tar_read(feature_partial_cor_ici_ratstamina)
  # partial_correlations = tar_read(feature_partial_cor_ici_nsclc)
  # annotations = tar_read(feature_annotations)
  # compound_type = "pathway"
  # compound_mapping = tar_read(metabolite_kegg)
  # lipid_type = "categories"
  if (is.na(partial_correlations$pcor[["partial_cor"]][1])) {
    network_partitioning = list(
      partitions = tibble::tibble(id = NA, ratio = NA, n_features = NA),
      q_value = NA
    )
  } else {
    network_correlations = partial_correlations$pcor |>
      dplyr::filter(significant) |>
      dplyr::transmute(start_node = s1, end_node = s2, weight = partial_cor)

    if (grepl("yeast", partial_correlations$data_id)) {
      use_annotation = annotations$yeast@annotation_features
    } else if (grepl("rat", partial_correlations$data_id)) {
      network_correlations2 = dplyr::left_join(
        network_correlations,
        compound_mapping[, c("feature_id", "kegg_id2")],
        by = c("start_node" = "feature_id"),
        relationship = "many-to-many"
      ) |>
        dplyr::mutate(s_kegg = kegg_id2, kegg_id2 = NULL)
      network_correlations2 = dplyr::left_join(
        network_correlations2,
        compound_mapping[, c("feature_id", "kegg_id2")],
        by = c("end_node" = "feature_id"),
        relationship = "many-to-many"
      ) |>
        dplyr::mutate(e_kegg = kegg_id2)
      network_correlations = network_correlations2 |>
        dplyr::transmute(
          start_node = s_kegg,
          end_node = e_kegg,
          weight = weight
        )

      if (compound_type %in% "pathway") {
        use_annotation = annotations$kegg_pathway@annotation_features
      }
    } else if (grepl("adenocarcinoma", partial_correlations$data_id)) {
      use_annotation = annotations$human@annotation_features
    } else if (grepl("egfrgenotype", partial_correlations$data_id)) {
      use_annotation = annotations$mouse@annotation_features
    } else if (grepl("nsclc", partial_correlations$data_id)) {
      if (lipid_type %in% "classes") {
        use_annotation = annotations$nsclc_classes@annotation_features
      } else {
        use_annotation = annotations$nsclc_categories@annotation_features
      }
      average_correlations = average_emf_correlations(
        network_correlations,
        annotations$nsclc_emfs
      )
      network_correlations = average_correlations
    }

    network_partitioning = calculate_qratio(
      network_correlations,
      use_annotation,
      use_weights = use_weights
    )
  }

  network_partitioning$data_id = partial_correlations$data_id
  network_partitioning$method_id = partial_correlations$method_id
  network_partitioning$full_id = partial_correlations$full_id
  return(network_partitioning)
}

pcor_pvalue_z = function(pcor_values) {
  mean_pcor = mean(pcor_values)
  sd_pcor = sd(pcor_values)
  z_scores = (pcor_values - mean_pcor) / sd_pcor

  p_value = 2 * pnorm(-abs(z_scores))
  p_adjust = p.adjust(p_value, method = "BH")
  tibble::tibble(pvalue = p_value, padjust = p_adjust)
}

pcor_pvalue_extreme = function(pcor_value_df, p_cut = 0.05) {
  # pcor_value_df = pcor_long
  # p_cut = 0.05
  mean_partial = mean(pcor_value_df$partial_cor)
  pcor_value_df = pcor_value_df |>
    dplyr::arrange(partial_cor) |>
    dplyr::mutate(
      rank = rank(partial_cor),
      fraction = rank / max(rank),
      p_value = dplyr::case_when(
        partial_cor < mean_partial ~ fraction,
        partial_cor > mean_partial ~ 1 - fraction
      ),
      significant = p_value <= (p_cut / 2)
    )
  pcor_value_df
}

n_extreme = function(in_value_df) {
  in_value_df = in_value_df |>
    dplyr::mutate(
      rank = rank(partial_cor),
      p_value = 1 - ((max(rank) - rank) / max(rank))
    )
}


compare_significant_partial_cor = function(feature_significant_combined) {
  # tar_load(feature_significant_combined)
  split_dataset = split(
    feature_significant_combined,
    feature_significant_combined$data_id
  )

  compared_indices = function(n1, n2) {
    l_n1 = length(n1)
    l_n2 = length(n2)
    len_union = length(base::union(n1, n2))
    len_intersect = length(base::intersect(n1, n2))
    min_len = min(c(l_n1, l_n2))

    jaccard = len_intersect / len_union
    overlap = len_intersect / min_len
    min_union = min_len + len_union
    log_combined = (log(len_intersect) + log(min_union)) -
      (log(min_len) + log(len_union)) -
      log(2)
    combined = exp(log_combined)
    #combined = (len_intersect * (min_len + len_union)) / (min_len * len_union) / 2

    data.frame(
      n1 = l_n1,
      n2 = l_n2,
      min = min_len,
      jaccard = jaccard,
      overlap = overlap,
      combined = combined
    )
  }

  out_compared = purrr::map(split_dataset, \(in_data) {
    # in_data = split_dataset[[1]]
    split_method = split(in_data, in_data$method_id)
    ici_method = split_method[["ici"]]

    compare_ici = purrr::map(split_method, \(x) {
      tmp_compare = compared_indices(ici_method$start_end, x$start_end)
      tmp_compare$d1 = "ici"
      tmp_compare$d2 = x$method_id[1]
      tmp_compare
    }) |>
      purrr::list_rbind()

    compare_ici$data_id = in_data$data_id[1]
    compare_ici
  }) |>
    purrr::list_rbind()
  out_compared
}

get_significant_partial_cor = function(feature_partial_cor) {
  # feature_partial_cor = tar_read(feature_partial_cor_ici_yeast)
  if (is.na(feature_partial_cor$pcor[["partial_cor"]][1])) {
    return(NULL)
  } else {
    network_correlations = feature_partial_cor$pcor |>
      dplyr::filter(significant) |>
      dplyr::transmute(start_node = s1, end_node = s2, weight = partial_cor) |>
      dplyr::filter(weight > 0) |>
      dplyr::mutate(
        start_end = paste0(start_node, ":", end_node),
        method_id = feature_partial_cor$method_id,
        data_id = feature_partial_cor$data_id
      )
    return(network_correlations)
  }
}

calculate_feature_partial_cor_pvalues = function(feature_correlation) {
  # feature_correlation = tar_read(feature_correlation_AN002783)
  if (is.null(feature_correlation)) {
    return(NULL)
  }

  partial_cor_pvalues = purrr::map(
    feature_correlation$cor_vals,
    calculate_partial_cor_pvalues
  )
  feature_correlation$partial_cor = partial_cor_pvalues
  feature_correlation
}

calculate_partial_cor_pvalues = function(feature_cor) {
  # feature_data = tar_read(feature_correlation_ici_yeast)
  # feature_data = tar_read(feature_correlation_pearson_base_nozero_ratstamina)
  # feature_data = tar_read(feature_correlation_ici_completeness_ratstamina)
  # feature_cor = tar_read(feature_correlation_AN002783)$cor_vals[[1]]
  if (is.null(feature_cor)) {
    return(NULL)
  }

  use_cor = feature_cor
  diag(use_cor) = 1
  use_cor[is.na(use_cor)] = 0
  pcor_vals = try(cor_to_pcor(use_cor))

  if (inherits(pcor_vals, "try-error")) {
    return(NULL)
  }
  pcor_long = cor_matrix_2_long_df(pcor_vals)
  pcor_long = pcor_long |>
    dplyr::filter(!is.na(cor)) |>
    dplyr::mutate(partial_cor = cor, cor = NULL, id = paste0(s1, ".", s2))

  feature_cor_df = cor_matrix_2_long_df(feature_cor) |>
    dplyr::mutate(id = paste0(s1, ".", s2))
  pcor_long = dplyr::left_join(
    feature_cor_df[, c("cor", "id")],
    pcor_long,
    by = "id"
  )
  pcor_long = pcor_long |>
    dplyr::filter(!(s1 == s2))
  pcor_long = pcor_long |>
    dplyr::filter(!is.na(partial_cor))
  pcor_p_values = pcor_pvalue_extreme(pcor_long)

  return(pcor_p_values)
  #pcor_vals[upper.tri(pcor_vals)] = NA
}

calculate_feature_correlation = function(
  keep_smd,
  predicted_annotation_datasets,
  n_bioreps = 20
) {
  # keep_smd = tar_read(smd_AN002783)
  # n_bioreps = 20

  if (is.null(keep_smd)) {
    return(NULL)
  }

  smd_metadata = metadata(keep_smd)
  if (!(smd_metadata$CHECK$ID %in% predicted_annotation_datasets$id)) {
    return(NULL)
  }

  if ((ncol(keep_smd) < n_bioreps)) {
    return(NULL)
  }

  feature_correlations = run_cor_everyway_new(keep_smd, features = TRUE)
  return(feature_correlations)
}

create_correlation_networks = function(
  counts_info,
  keep_num,
  sample_col,
  class_col
) {
  # counts_info = tar_read(yeast_counts_info)
  # keep_num = 1
  # sample_col = "sample"
  # class_col = "treatment"
  counts = counts_info$counts
  info = counts_info$info

  if (length(class_col) == 2) {
    filter_col = class_col[1]
    median_col = class_col[2]
  } else {
    filter_col = class_col
    median_col = class_col
  }
  counts_filter = t(keep_non_zero_percentage(
    t(counts),
    sample_classes = info[[filter_col]],
    keep_num = keep_num
  ))
  counts_completeness = pairwise_completeness(counts_filter)
  counts_cor = run_cor_everyway(t(counts_filter), counts_completeness)
  counts_cor
}


compound_annotation = function(data_file, type = "pathway") {
  # data_file = "data/kegg_compound_mapping.rds"
  # type = "pathway"
  compound_data = readRDS(data_file)

  use_entries = compound_data[[type]]
  split_data = split(use_entries$compound, use_entries$id) |>
    purrr::map(unique)

  pathway_meta = compound_data$meta |>
    dplyr::filter(type %in% type)
  pathway_description = pathway_meta[["description"]]
  names(pathway_description) = pathway_meta$id
  pathway_description = pathway_description[names(split_data)]
  annotation_obj = categoryCompare2::annotation(
    split_data,
    annotation_type = paste0("kegg-", type),
    description = pathway_description,
    feature_type = "kegg-compound"
  )
  annotation_obj
}

get_feature_annotations = function(kegg_data, nsclc_lipidclasses) {
  #tar_load(kegg_data)
  #tar_load(nsclc_lipidclasses)
  compound_pathway = compound_annotation(kegg_data, "pathway")
  compound_network = compound_annotation(kegg_data, "network")
  compound_module = compound_annotation(kegg_data, "module")

  yeast_keys = AnnotationDbi::keys(org.Sc.sgd.db)
  yeast_entrez = suppressMessages(AnnotationDbi::select(
    org.Sc.sgd.db,
    keys = yeast_keys,
    columns = "ENTREZID"
  ))
  yeast_uniq_entrez = unique(yeast_entrez$ENTREZID)
  yeast_uniq_entrez = yeast_uniq_entrez[!is.na(yeast_uniq_entrez)]
  yeast_reactome = suppressMessages(AnnotationDbi::select(
    reactome.db,
    keys = yeast_uniq_entrez,
    keytype = "ENTREZID",
    columns = c("PATHNAME", "REACTOMEID")
  ))
  yeast_entrez = yeast_entrez |>
    dplyr::transmute(ID = ORF, ENTREZID = ENTREZID)
  yeast_id_reactome = dplyr::left_join(
    yeast_entrez,
    yeast_reactome,
    by = "ENTREZID"
  )

  yeast_annotation = create_feature_annotation_object(yeast_id_reactome)

  mouse_keys = AnnotationDbi::keys(org.Mm.eg.db)
  mouse_entrez = suppressMessages(AnnotationDbi::select(
    org.Mm.eg.db,
    keys = mouse_keys,
    columns = "ENSEMBL"
  ))
  mouse_uniq_entrez = unique(mouse_entrez$ENTREZID)
  mouse_uniq_entrez = mouse_uniq_entrez[!is.na(mouse_uniq_entrez)]
  mouse_reactome = suppressMessages(AnnotationDbi::select(
    reactome.db,
    keys = mouse_uniq_entrez,
    keytype = "ENTREZID",
    columns = c("PATHNAME", "REACTOMEID")
  ))
  mouse_entrez = mouse_entrez |>
    dplyr::transmute(ID = ENSEMBL, ENTREZID = ENTREZID)
  mouse_id_reactome = dplyr::left_join(
    mouse_entrez,
    mouse_reactome,
    by = "ENTREZID",
    relationship = "many-to-many"
  )
  mouse_annotation = create_feature_annotation_object(mouse_id_reactome)

  human_keys = AnnotationDbi::keys(org.Hs.eg.db)
  human_entrez = suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db,
    keys = human_keys,
    columns = "ENSEMBL"
  ))
  human_uniq_entrez = unique(human_entrez$ENTREZID)
  human_uniq_entrez = human_uniq_entrez[!is.na(human_uniq_entrez)]
  human_reactome = suppressMessages(AnnotationDbi::select(
    reactome.db,
    keys = human_uniq_entrez,
    keytype = "ENTREZID",
    columns = c("PATHNAME", "REACTOMEID")
  ))
  human_entrez = human_entrez |>
    dplyr::transmute(ID = ENSEMBL, ENTREZID = ENTREZID)
  human_id_reactome = dplyr::left_join(
    human_entrez,
    human_reactome,
    by = "ENTREZID",
    relationship = "many-to-many"
  )
  human_annotation = create_feature_annotation_object(human_id_reactome)

  nsclc_annotation_df = readRDS(nsclc_lipidclasses)
  nsclc_annotation_df_categories = nsclc_annotation_df |>
    dplyr::mutate(ID = emf, REACTOMEID = Voted.Categories)
  nsclc_annotation_categories = create_feature_annotation_object(
    nsclc_annotation_df_categories
  )

  nsclc_annotation_df_classes = nsclc_annotation_df |>
    dplyr::mutate(ID = emf, REACTOMEID = Voted.Classes)
  nsclc_annotation_classes = create_feature_annotation_object(
    nsclc_annotation_df_classes
  )

  nsclc_emfs = nsclc_annotation_df |>
    dplyr::select(feature_id, sudo_EMF, emf) |>
    dplyr::distinct()

  list(
    human = human_annotation,
    mouse = mouse_annotation,
    yeast = yeast_annotation,
    nsclc_categories = nsclc_annotation_categories,
    nsclc_classes = nsclc_annotation_classes,
    nsclc_emfs = nsclc_emfs,
    kegg_pathway = compound_pathway,
    kegg_module = compound_module,
    kegg_network = compound_network
  )
}

create_feature_annotation_object = function(annotation_df) {
  # annotation_df = yeast_id_reactome
  annotation_df = annotation_df |>
    dplyr::filter(!is.na(REACTOMEID))
  split_data = split(annotation_df$ID, annotation_df$REACTOMEID) |>
    purrr::map(unique)

  if (!is.null(annotation_df[["PATHNAME"]])) {
    pathway_meta = annotation_df |>
      dplyr::select(PATHNAME, REACTOMEID) |>
      dplyr::distinct()
    pathway_description = pathway_meta[["PATHNAME"]]
    names(pathway_description) = pathway_meta[["REACTOMEID"]]
    pathway_description = pathway_description[names(split_data)]
    annotation_obj = categoryCompare2::annotation(
      split_data,
      description = pathway_description
    )
  } else {
    annotation_obj = categoryCompare2::annotation(split_data)
  }

  annotation_obj
}

average_emf_correlations = function(nsclc_partial_cor, nsclc_emf_mapping) {
  # nsclc_partial_cor = network_correlations
  # nsclc_emf_mapping = annotations$nsclc_emfs
  nsclc_partial_cor = nsclc_partial_cor |>
    dplyr::filter(
      (start_node %in% nsclc_emf_mapping$feature_id) &
        (end_node %in% nsclc_emf_mapping$feature_id)
    )
  nsclc_partial_cor_emf = dplyr::left_join(
    nsclc_partial_cor,
    nsclc_emf_mapping |>
      dplyr::transmute(feature_id = feature_id, start_emf = emf),
    by = c("start_node" = "feature_id"),
    relationship = "many-to-many"
  )
  nsclc_partial_cor_emf = dplyr::left_join(
    nsclc_partial_cor_emf,
    nsclc_emf_mapping |>
      dplyr::transmute(feature_id = feature_id, end_emf = emf),
    by = c("end_node" = "feature_id"),
    relationship = "many-to-many"
  )
  cross_emf = nsclc_partial_cor_emf |>
    dplyr::filter(!(start_emf == end_emf)) |>
    dplyr::mutate(cross_emf_id = paste0(start_emf, ".", end_emf))
  average_cross_emf = cross_emf |>
    dplyr::group_by(cross_emf_id) |>
    dplyr::summarise(weight2 = sqrt(mean(weight^2))) |>
    dplyr::mutate(weight = weight2)
  split_emfs = stringr::str_split_fixed(
    average_cross_emf$cross_emf_id,
    "\\.",
    2
  )
  average_cross_emf$start_node = split_emfs[, 1]
  average_cross_emf$end_node = split_emfs[, 2]
  average_cross_emf = average_cross_emf |>
    dplyr::select(start_node, end_node, weight)
  average_cross_emf
}

cleanup_qratio_table = function(feature_qratio_summary) {
  # tar_load(feature_qratio_summary)
  map_method = tibble::tribble(
    ~method_id            , ~Method ,
    "ici"                 , "IK"    ,
    "ici_completeness"    , "IKC"   ,
    "pearson_base"        , "PB"    ,
    "pearson_base_nozero" , "PN0"   ,
    "pearson_log1p"       , "PL1"   ,
    "pearson_log"         , "PL"    ,
    "kt"                  , "Kt"
  )
  map_dataset = tibble::tribble(
    ~data_id                      , ~Dataset                   ,
    "barton_yeast"                , "Yeast RNA-Seq"            ,
    "brainsonrnaseq_egfrgenotype" , "EGFR Genotype RNA-Seq"    ,
    "mwtab_ratstamina"            , "Rat Stamina Metabolomics" ,
    "nsclc"                       , "NSCLC Lipidomics"
  )
  feature_qratio_summary = dplyr::left_join(
    feature_qratio_summary,
    map_method,
    by = "method_id"
  )
  feature_qratio_summary = dplyr::left_join(
    feature_qratio_summary,
    map_dataset,
    by = "data_id"
  )
  feature_qratio_out = feature_qratio_summary |>
    dplyr::ungroup() |>
    dplyr::transmute(Dataset, Method, Q = q_value)
  feature_qratio_wide = feature_qratio_out |>
    tidyr::pivot_wider(
      id_cols = Method,
      values_from = "Q",
      names_from = Dataset
    )
  feature_qratio_ft = flextable::flextable(feature_qratio_wide)
  feature_qratio_ft = flextable::colformat_double(feature_qratio_ft, digits = 3)

  highlight_locs = matrix(c(1, 2, 1, 3, 2, 4, 2, 5), nrow = 4, byrow = TRUE)
  feature_qratio_highlight = feature_qratio_ft
  for (irow in seq_len(nrow(highlight_locs))) {
    feature_qratio_highlight = flextable::bold(
      feature_qratio_highlight,
      i = highlight_locs[irow, 1],
      j = highlight_locs[irow, 2],
      part = "body"
    )
  }
  feature_qratio_highlight
}

cleanup_significant_compare_table = function(in_summary) {
  # in_summary = tar_read(feature_compare_summary)
  in_summary = in_summary |>
    dplyr::filter(!(d1 == d2))
  flex_summary = in_summary |>
    flextable::flextable() |>
    flextable::colformat_double(j = c(4, 5, 6), digits = 3) |>
    flextable::set_table_properties(layout = "autofit")
  flex_summary
}
