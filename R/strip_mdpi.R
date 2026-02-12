strip_headers = function(
  in_rmd = "docs/ici_kt_manuscript.Rmd",
  out_rmd = "docs/ici_kt_manuscript_mdpi.Rmd"
) {
  # in_rmd = "docs/ici_kt_manuscript.Rmd"
  # out_rmd = "docs/ici_kt_manuscript_mdpi.Rmd"

  rmd_doc = readLines(in_rmd)

  has_header = grepl("^##+ ", rmd_doc)
  has_bullet = grepl("^\\* ", rmd_doc)
  doc_out = rmd_doc
  doc_out[has_header] = gsub("^##+ ", "", rmd_doc[has_header])
  doc_out[has_bullet] = gsub("^\\* ", "", rmd_doc[has_bullet])
  cat(doc_out, file = out_rmd, sep = "\n")
  out_rmd
}
