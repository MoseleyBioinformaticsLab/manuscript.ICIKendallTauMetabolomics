args = commandArgs(trailingOnly = TRUE)
# print(args[1])
# args = c(1)
options(future.rng.onMisuse = "ignore")
library(furrr)
plan(multicore, workers = as.numeric(args[1]))
library(ICIKendallTau)
out_log_file = paste0("python/r_log_", args[1], ".log")
enable_logging(out_log_file, memory = TRUE)

n_feature = 10000
n_sample = 400

set.seed(1234)
in_data = rnorm(n_feature * n_sample) |>
  matrix(nrow = n_feature, ncol = n_sample)
colnames(in_data) = paste0("s", seq_len(n_sample))
ICIKendallTau:::log_memory()

ici_res = ici_kendalltau(in_data, return_matrix = FALSE)

out_csv_file = paste0("python/r_res_", args[1], ".csv")
write.table(
  c(args[1], ici_res$run_time, out_log_file),
  file = out_csv_file,
  row.names = FALSE,
  col.names = FALSE,
  sep = ","
)
