fgsea_run <- function(vd, term2feature, term_info, min_size = 3) {
  vd <- vd %>%
    filter(!is.na(value) & !is.na(feature_id))
  ranks <-  set_names(vd$value, vd$feature_id)
  fgsea::fgsea(pathways = term2feature, stats = ranks, nproc = 6, minSize = min_size, eps = 0) %>%
    as_tibble() %>%
    left_join(term_info, by = c("pathway" = "term_id")) %>%
    arrange(NES) %>%
    select(term_id = pathway, term_name, pval, padj, NES, size, leading_edge = leadingEdge)
}



fgsea_all_terms <- function(d, all_terms, feature_var = "gene_symbol", value_var = "logFC", group_var = "contrast") {
  val_data <- d |>
    rename(value = !!value_var, feature_id = !!feature_var, group = !!group_var)
  ontologies <- names(all_terms)
  map_dfr(ontologies, function(ont) {
    cat(str_glue("  Computing fgsea for {ont}\n\n"))
    terms <- all_terms[[ont]]
    term2feature <- terms$mapping |>
      rename(feature_id = !!feature_var) |>
      group_by(term_id) |>
      summarise(features = list(feature_id)) |>
      deframe()
    val_data |>
      group_split(group) |>
      map_dfr(function(w) {
        fgsea_run(w, term2feature, terms$terms) |>
          add_column(!!group_var := first(w$group))
      }) |>
      add_column(ontology = ont)
  })
}


tabulate_gse <- function(gse, fdr_limit = 1) {
  gse |>
    filter(padj < fdr_limit) |>
    unnest(leading_edge) |>
    group_by(term_id, term_name, NES, pval, padj, ontology, contrast) |>
    summarise(genes = str_c(leading_edge, collapse = ", ")) |>
    arrange(NES) |>
    ungroup() |>
    relocate(ontology, .before = 1)
}



plot_fgsea_enrichment <- function(term, res, term_data, value = "logFC", feature = "gene_name") {
  lst <- term_data$term2feature[[term]]
  rnks <- set_names(res[[value]], res[[feature]])
  fgsea::plotEnrichment(lst, rnks)
}

