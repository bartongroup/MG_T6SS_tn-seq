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
  add_row(colname = "PFAMs", prefix = "") |>
  add_column(ontology =  c(rep("KEGG", nrow(KEGG_COLUMNS)), "GO", "PFAM"))


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


get_go_terms <- function(obo_file = "http://purl.obolibrary.org/obo/go.obo") {
  parsed <- readr::read_lines(obo_file) |>
    parse_obo_file()

  terms <- parsed |>
    dplyr::filter(key == "name") |>
    select(term_id, term_name = value)

  alt_terms <- parsed |>
    dplyr::filter(key == "alt_id") |>
    left_join(terms, by = "term_id") |>
    select(term_id = value, term_name)

  bind_rows(
    terms,
    alt_terms
  )
}


parse_obo_file <- function(obo) {
  # Index start and end of each term
  idx_start_term <- obo |>
    stringr::str_which("\\[Term\\]")
  idx_empty <- obo |>
    stringr::str_which("^$")
  idx_empty <- idx_empty[idx_empty > idx_start_term[1]]
  idx_end_term <- idx_empty[1:length(idx_start_term)]

  purrr::map2_dfr(idx_start_term, idx_end_term, function(i1, i2) {
    obo_term <- obo[(i1 + 1):(i2 - 1)]

    trm <- obo_term |>
      stringr::str_split(":\\s", 2, simplify = TRUE)
    colnames(trm) <- c("key", "value")
    # assuming term_id is in the first line
    tid <- trm[1, 2]
    cbind(trm, term_id = tid) |>
      as.data.frame()
  }) |>
    tibble::as_tibble()
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
  # I discovered that pathways koxxxx and mapxxx are identical, so removing all
  # mapxxx terms.

  egmap <- egmap |>
    filter(!(str_detect(term_id, "^map") & colname == "KEGG_Pathway"))

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


eggnog_pfam <- function(egmap) {
  mapping <- egmap |>
    filter(ontology == "PFAM") |>
    select(feature_id = id, term_id) |>
    distinct()
  terms <- tibble(
    term_id = unique(mapping$term_id),
    term_name = term_id
  )

  list(
    terms = terms,
    mapping = mapping
  )
}


