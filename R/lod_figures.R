graph_median_min = function(median_correlation)
{
  # median_correlation = tar_read(correlate_ranks_yeast)
  median_plots = create_median_plot(median_correlation$medians, median_correlation$correlation)
  min_plots = create_min_median_plot(median_correlation$medians, median_correlation$correlation)
  list(median = median_plots,
       min = min_plots)
}

find_good_location = function(all_values, perc_range = 0.1)
{
  # range_values = c(-0.8652091, 18.1403589)
  # perc_range = 0.1
  
  range_values = range(all_values, na.rm = TRUE)
  range_diff = range_values[2] - range_values[1]
  
  use_loc = range_values[1] + (range_diff * perc_range)
  use_loc
}

calculate_median_min_correlation = function(na_ranks)
{
  # na_ranks = tar_read(rank_ordered_yeast)
  med_ranks = na_ranks$ranks |>
    dplyr::group_by(n_na, treatment) |>
    dplyr::summarise(median_median = median(median_rank),
                     median_min = min(median_rank)) |>
    dplyr::ungroup()
  
  med_cor = med_ranks|>
    dplyr::group_by(treatment) |>
    dplyr::summarise(kt_median = ici_kt(n_na, median_median)[1],
                     kt_min = ici_kt(n_na, median_min)[1])
  return(list(medians = med_ranks,
              correlation = med_cor,
            id = na_ranks$id))
  
}

create_median_plot = function(median_stuff, correlation_stuff)
{
  # tmp = tar_read(correlate_ranks_yeast)
  # median_stuff = tmp$medians
  # correlation_stuff = tmp$correlation
  n_treatments = length(unique(median_stuff$treatment))
  if (n_treatments <= 3) {
    x_loc_fraction = 0.7
  } else {
    x_loc_fraction = 0.6
  }
  use_locations = median_stuff |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(cor_x = find_good_location(median_median, x_loc_fraction),
                     cor_y = find_good_location(n_na, 0.95))
  use_locations = dplyr::left_join(use_locations, correlation_stuff, by = "treatment")
  use_locations = use_locations |>
    dplyr::mutate(cor_label = paste0("\u03C4: ", format(kt_median, digits = 2)))
  
  nrow = 1
  if (n_treatments > 5) {
    nrow = 3
  } else if (n_treatments > 3) {
    nrow = 2
  }
  out_plot = median_stuff |>
    ggplot(aes(x = median_median, y = n_na)) +
    geom_point() +
    geom_label(data = use_locations, aes(x = cor_x, y = cor_y, label = cor_label), hjust = 0) +
    facet_wrap(~ treatment, nrow = nrow, scales = "free") +
    theme(strip.background = NULL,
          strip.text.x = element_text(hjust = 0)) +
    labs(x = "Median(Median Rank)", y = "N-Missing")
  
  out_plot
}

create_min_median_plot = function(median_stuff, correlation_stuff)
{
  # tmp = tar_read(correlate_ranks_yeast)
  # median_stuff = tmp$medians
  # correlation_stuff = tmp$correlation
  # 
  # tmp = tar_read(correlate_medians_typeandtumorculture)
  # median_stuff = tmp$median_min
  # xlimit = tmp$quantile
  n_treatments = length(unique(median_stuff$treatment))
  if (n_treatments <= 3) {
    x_loc_fraction = 0.7
  } else {
    x_loc_fraction = 0.6
  }
  
  use_locations = median_stuff |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(cor_x = find_good_location(median_min, x_loc_fraction),
                     cor_y = find_good_location(n_na, 0.95))
  use_locations = dplyr::left_join(use_locations, correlation_stuff, by = "treatment")
  use_locations = use_locations |>
    dplyr::mutate(cor_label = paste0("\u03C4: ", format(kt_min, digits = 2)))
  
  nrow = 1
  if (n_treatments > 5) {
    nrow = 3
  } else if (n_treatments > 3) {
    nrow = 2
  }
  
  out_plot = median_stuff |>
    ggplot(aes(x = median_min, y = n_na)) +
    geom_point() +
    geom_label(data = use_locations, aes(x = cor_x, y = cor_y, label = cor_label), hjust = 0) +
    facet_wrap(~ treatment, nrow = nrow, scales = "free") +
    theme(strip.background = NULL,
          strip.text.x = element_text(hjust = 0)) +
    labs(x = "Min(Median-Rank)", y = "N-Missing")
  
  out_plot
}
