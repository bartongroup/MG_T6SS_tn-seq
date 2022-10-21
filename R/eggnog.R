KEGG_BASE_URL <- "https://rest.kegg.jp"

KEGG_COLUMNS <- tribble(
  ~colname, ~dbentry, ~prefix,
  "KEGG_ko", "ko", "",
  "KEGG_Pathway", "pathway/ko", "",
  "KEGG_Pathway", "pathway/map", "",
  "KEGG_Module", "module", "md:",
  "KEGG_Reaction", "reaction", "rn:",
  "KEGG_rclass", "rclass", "rc:",
  "BRITE", "brite", "br:"
)

EGGNOG_COLUMNS <- KEGG_COLUMNS |>
  add_row(colname = "GOs", prefix = "") |>
  add_column(ontology =  c(rep("KEGG", nrow(KEGG_COLUMNS)), "GO"))


fetch_kegg_data <- function(path) {
  url <- file.path(KEGG_BASE_URL, path)
  res <- httr::GET(url)
  rawToChar(res$content)
}

fetch_kegg_pathways <- function(dbentries) {
  query <- stringr::str_glue("list/{dbentries}")
  s <- fetch_kegg_data(query)
  readr::read_tsv(I(s), col_names = c("term_id", "term_name"), show_col_types = FALSE) |>
    dplyr::mutate(term_id = stringr::str_remove(term_id, "path:"))
}

get_kegg_terms <- function() {
  KEGG_COLUMNS |>
    rowwise() |>
    group_split() |>
    map_dfr(function(r) {
      fetch_kegg_pathways(r$dbentry) |>
        add_column(colname = r$colname)
    })
}

extract_obo_values <- function(obo, key) {
  obo |>
    stringr::str_subset(stringr::str_glue("^{key}:")) |>
    stringr::str_remove(stringr::str_glue("{key}:\\s"))
}

get_go_terms <- function(obo_file = "http://purl.obolibrary.org/obo/go.obo") {
  obo <- readr::read_lines(obo_file)
  ids <- extract_obo_values(obo, "id")
  names <- extract_obo_values(obo, "name")

  tibble::tibble(
    term_id = ids,
    term_name = names
  )
}


##############################################################

read_eggnog_annotation_file <- function(file) {
  read_tsv(file, comment = "##", show_col_types = FALSE) |>
    rename(id = `#query`)
}


read_eggnog_mapping <- function(file) {
  raw <- read_eggnog_annotation_file(file)
  EGGNOG_COLUMNS |>
    rowwise() |>
    group_split() |>
    map_dfr(function(r) {
      raw[, c("id", r$colname)] |>
        rename(term_id = !!r$colname) |>
        filter(term_id != "-") |>
        separate_rows(term_id, sep = ",") |>
        mutate(term_id = paste0(r$prefix, term_id)) |>
        add_column(colname = r$colname, ontology = r$ontology)
    })
}

eggnog_kegg <- function(egmap) {
  kg <- get_kegg_terms()

  colnames <- egmap |>
    filter(ontology == "KEGG") |>
    pull(colname) |>
    unique()

  map(colnames, function(name) {
    mapping <- egmap |>
      filter(colname == name) |>
      select(feature_id = id, term_id) |>
      distinct()
    terms <- kg |>
      filter(colname == name)
    list(
      terms = terms,
      mapping = mapping
    )
  }) |>
    set_names(colnames)
}


eggnog_go <- function(egmap) {
  go <- get_go_terms()

  mapping <- egmap |>
    filter(ontology == "GO") |>
    select(feature_id = id, term_id) |>
    distinct()

  list(
    terms = go,
    mapping = mapping
  )
}


