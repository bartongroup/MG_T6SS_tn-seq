prepare_terms <- function(all_terms) {
  ontologies <- names(all_terms)
  map(ontologies, function(ont) {
    tms <- all_terms[[ont]]
    fenr::prepare_for_enrichment(tms$terms, tms$mapping, feature_name = "feature_id")
  }) |>
    set_names(ontologies)
}


write_rds_name <- function(obj) {
  path <- file.path("shiny", "data")
  if (!dir.exists(path)) dir.create(path)
  obj_name <- deparse(substitute(obj))
  file_name <- file.path(path, str_glue("{obj_name}.rds"))
  write_rds(obj, file_name, compress = "xz")
}

save_data_for_shiny <- function(cnts, edger, fterms, gse) {
  write_rds_name(cnts)
  write_rds_name(edger)
  write_rds_name(fterms)
  write_rds_name(gse)
}

