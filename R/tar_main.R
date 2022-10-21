targets_main <- function() {


  read_data <- list(
    tar_target(metadata, make_metadata()),
    tar_target(cnts, read_counts(COUNT_FILE, metadata, MIN_COUNT))
  )

  differential <- list(
    tar_target(edger, edger_local(cnts, CONTRASTS)),
    tar_target(gse, fgsea_all_terms(edger, all_terms, feature_var = "feature_id"))
  )

  eggnog <- list(
    tar_target(egg, read_eggnog_mapping(EGGNOG_FILE)),
    tar_target(egg_kegg, eggnog_kegg(egg)),
    tar_target(egg_go, eggnog_go(egg)),
    tar_target(all_terms, c(list(GO = egg_go), egg_kegg))
  )

  report <- list(
    tar_target(min_count, MIN_COUNT),
    tar_target(fig_count_dist, plot_count_distribution(cnts, MIN_COUNT)),
    tar_target(fig_replicates, plot_replicate_comparison(cnts)),
    tar_target(fig_volcano, plot_volcano(edger)),
    tar_target(fig_ma, plot_ma(edger)),
    tar_target(tab_gse, tabulate_gse(gse)),
    tar_target(example_genes, get_extreme_genes(edger, 1, 8)),
    tar_target(fig_example_genes, plot_genes(cnts, example_genes))
  )

  info <- list(
    tar_target(edger_version, packageVersion("edgeR")),
    tar_target(fgsea_version, packageVersion("fgsea"))
  )

  c(
    read_data,
    differential,
    eggnog,
    info,
    report
  )

}
