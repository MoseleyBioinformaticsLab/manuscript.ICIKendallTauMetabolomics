compare_previous_new_samples = function() {
  old_dir = "/big_data/data/rmflight_ici-kendallt/_targets/objects"
  new_dir = "/big_data/data/rmflight_icikt_metabolomics/_targets/objects"

  old_env = new.env()
  new_env = new.env()

  old_env$realistic_sample_1 = readRDS(fs::path(old_dir, "realistic_sample_1"))
  old_env$realistic_sample_2 = readRDS(fs::path(old_dir, "realistic_sample_2"))

  old_data = tibble::tibble(
    s1 = old_env$realistic_sample_1,
    s2 = old_env$realistic_sample_2
  )

  new_env$realistic_sample_1 = readRDS(fs::path(new_dir, "realistic_sample_1"))
  new_env$realistic_sample_2 = readRDS(fs::path(new_dir, "realistic_sample_2"))

  new_data = tibble::tibble(
    s1 = new_env$realistic_sample_1,
    s2 = new_env$realistic_sample_2
  )
}
