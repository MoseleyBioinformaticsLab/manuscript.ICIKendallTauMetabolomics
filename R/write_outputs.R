write_missingness_table = function(
  missingness_test_table,
  file = "docs/supplementary_missingness_table.xlsx"
) {
  openxlsx::write.xlsx(missingness_test_table, file = file)
}
