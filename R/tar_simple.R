targets_simple <- function() {

  read_data <- list(
    tar_target(metadata, make_metadata()),
    tar_target(cnts, read_counts(COUNT_FILE, metadata)),
    tar_target(n_genes, length(cnts$sel))
  )

  differential_expression <- list(
    tar_target(edger, edger_local(cnts, CONTRASTS, logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT)),

    tar_target(fig_volcano, plot_volcano(edger)),
    tar_target(fig_ma, plot_ma(edger)),

    tar_target(sig_genes, edger |> filter(sig) |> pull(feature_id) |> unique()),
    tar_target(fig_sig_genes, plot_genes(cnts, sig_genes))
  )

  report <- list(
    tar_target(fig_count_dist, plot_count_distribution(cnts)),
    tar_target(fig_replicates, plot_replicate_comparison(cnts))
  )

  c(
    read_data,
    differential_expression,
    report
  )

}
