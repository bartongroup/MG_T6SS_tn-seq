targets_main <- function() {


  read_data <- list(
    tar_target(metadata, make_metadata()),
    tar_target(cnts, read_counts(COUNT_FILE, metadata))
  )

  eggnog <- list(
    tar_target(egg, read_eggnog_mapping(EGGNOG_FILE)),
    tar_target(egg_kegg, eggnog_kegg(egg)),
    tar_target(egg_go, eggnog_go(egg)),
    tar_target(all_terms, c(list(GO = egg_go), egg_kegg)),
    tar_target(fterms, prepare_terms(all_terms))
  )

  map_cv <- tar_map(
    values = CV_LIMITS,
    names = LIMIT,

    tar_target(cnts_flt, filter_cv(cnts, cv_limit = LIMIT)),
    tar_target(n_genes, length(cnts_flt$sel)),

    tar_target(edger, edger_local(cnts_flt, CONTRASTS, logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT)),
    tar_target(gse, fgsea_all_terms(edger, all_terms, feature_var = "feature_id")),

    tar_target(fig_volcano, plot_volcano(edger)),
    tar_target(fig_ma, plot_ma(edger)),
    tar_target(tab_gse, tabulate_gse(gse)),

    tar_target(sig_genes, edger |> filter(sig) |> pull(feature_id) |> unique()),
    tar_target(fig_sig_genes, plot_genes(cnts_flt, sig_genes))
  )

  report <- list(
    tar_target(fig_count_dist, plot_count_distribution(cnts)),
    tar_target(fig_replicates, plot_replicate_comparison(cnts)),
    tar_target(lysine_deg_genes, gse_0.5 |> filter(term_id == "ko00310" & contrast == "YL57-YL37") |> pull(leading_edge) |> unlist()),
    tar_target(fig_lysin_deg, plot_genes(cnts_flt_0.5, lysine_deg_genes)),
    tar_target(chaperones_genes, gse_0.5 |> filter(term_id == "br:ko03110" & contrast == "GM103-YL37") |> pull(leading_edge) |> unlist()),
    tar_target(fig_chaperones, plot_genes(cnts_flt_0.5, chaperones_genes)),

    tar_target(sav_de, write_tsv(edger_0.5, "tab/de_cv0.5.tsv")),
    tar_target(sav_gse, write_tsv(tab_gse_0.5, "tab/gse_cv0.5.tsv"))
  )

  shiny <- list(
    tar_target(sav_shiny, save_data_for_shiny(cnts_flt_0.5, edger_0.5, fterms, gse_0.5))
  )

  info <- list(
    tar_target(edger_version, packageVersion("edgeR")),
    tar_target(fgsea_version, packageVersion("fgsea"))
  )

  c(
    read_data,
    map_cv,
    eggnog,
    info,
    report,
    shiny
  )

}
