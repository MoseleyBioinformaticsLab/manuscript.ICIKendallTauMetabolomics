work_out_similarities = function() {
  test_annotation_data = read.table(
    "~/Projects/work/collaborators/nate_helsley/hepatocellularcarcinoma_rnaseq_metabolomics/data/ChEBI2Reactome_All_Levels.txt",
    sep = "\t"
  )
  test_annotation_data = test_annotation_data |>
    dplyr::filter(grepl("^R-HSA", V2))
  annotations_2_feature = split(
    test_annotation_data$V1,
    test_annotation_data$V2
  )

  annotations_2_feature = purrr::map(annotations_2_feature, unique)
  keep_annotation = purrr::map_lgl(annotations_2_feature, \(x) {
    (length(x) >= 15) & (length(x <= 500))
  })

  annotations_2_feature = annotations_2_feature[keep_annotation]
  n_do = length(annotations_2_feature)
  similarity_matrix = matrix(
    NA,
    n_do,
    n_do
  )

  for (ianno in seq(1, n_do)) {
    for (janno in seq(ianno, n_do)) {
      similarity_matrix[ianno, janno] = similarity_matrix[
        janno,
        ianno
      ] = categoryCompare2:::combined_coefficient(
        annotations_2_feature[[ianno]],
        annotations_2_feature[[janno]]
      )
    }
  }

  rownames(similarity_matrix) = colnames(similarity_matrix) = names(
    annotations_2_feature
  )

  high_similarity = similarity_matrix
  high_similarity[similarity_matrix < 0.6] = 0

  adj_graph = igraph::graph_from_adjacency_matrix(
    high_similarity,
    mode = "upper",
    weighted = TRUE,
    diag = FALSE
  )

  adj_communities = igraph::cluster_walktrap(adj_graph)
  adj_membership = igraph::membership(adj_communities)

  membership_list = split(names(adj_membership), adj_membership)

  collapsed_features = purrr::map(membership_list, function(x) {
    # x = membership_list[[1]]
    all_annote = annotations_2_feature[x] |>
      unlist(use.names = FALSE) |>
      unique()
  })

  n_collapsed = length(collapsed_features)
  collapsed_similarity = matrix(NA, n_collapsed, n_collapsed)

  for (icol in seq(1, n_collapsed)) {
    for (jcol in seq(icol, n_collapsed)) {
      collapsed_similarity[icol, jcol] = categoryCompare2::combined_coefficient(
        collapsed_features[[icol]],
        collapsed_features[[jcol]]
      )
    }
  }
  collapsed_similarity[lower.tri(collapsed_similarity, diag = TRUE)] = NA
}

determine_similar_annotations = function(all_feature_annotations) {
  split_feature_annotations = split(
    all_feature_annotations$features,
    all_feature_annotations$annotation
  )
  split_feature_annotations = purrr::map(split_feature_annotations, unique)

  keep_annotation = purrr::map_lgl(split_feature_annotations, \(x) {
    (length(x) >= 15) & (length(x) <= 500)
  })

  split_feature_annotations = split_feature_annotations[keep_annotation]
  n_do = length(split_feature_annotations)
  similarity_matrix = matrix(
    NA,
    n_do,
    n_do
  )

  for (ianno in seq(1, n_do)) {
    for (janno in seq(ianno, n_do)) {
      similarity_matrix[ianno, janno] = similarity_matrix[
        janno,
        ianno
      ] = categoryCompare2:::combined_coefficient(
        split_feature_annotations[[ianno]],
        split_feature_annotations[[janno]]
      )
    }
  }

  rownames(similarity_matrix) = colnames(similarity_matrix) = names(
    split_feature_annotations
  )

  high_similarity = similarity_matrix
  high_similarity[similarity_matrix < 0.6] = 0

  adj_graph = igraph::graph_from_adjacency_matrix(
    high_similarity,
    mode = "upper",
    weighted = TRUE,
    diag = FALSE
  )

  adj_communities = igraph::cluster_walktrap(adj_graph)
  adj_membership = igraph::membership(adj_communities)

  membership_list = split(names(adj_membership), adj_membership)
}


get_all_annotations = function(file_locs) {
  # file_locs = predicted_annotation_datasets$dir

  all_feature_annotations = purrr::map(
    file_locs,
    get_dataset_annotation,
    .progress = TRUE
  ) |>
    purrr::list_rbind()

  all_feature_annotations
}

get_dataset_annotation = function(file_dir) {
  # file_dir = "predicted_annotations/AN000001"

  json_loc = fs::path(file_dir, "xref", "reactome", "predictions-union.json")

  json_annotation = categoryCompare2::json_2_annotation(json_loc)
  annotation_features = json_annotation@annotation_features
  multi_species_single = tibble::tibble(species = names(annotation_features)) |>
    dplyr::mutate(nonspecies = stringr::str_replace_all(species, "-(.*)-", "-"))
  split_nonspecies = split(
    multi_species_single$species,
    multi_species_single$nonspecies
  )

  combined_annotation = purrr::map(split_nonspecies, \(in_split) {
    annotation_features[in_split] |> unlist(use.names = FALSE) |> unique()
  })

  annotation_tbl = purrr::imap(combined_annotation, \(features, annotation) {
    tibble::tibble(
      features = features,
      annotation = annotation,
      id = fs::path_file(file_dir)
    )
  }) |>
    purrr::list_rbind()

  return(annotation_tbl)
}
