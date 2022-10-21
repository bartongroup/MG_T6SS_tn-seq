edger_local <- function(dset, contrasts, fdr_limit = 0.05, logfc_limit = 0) {
  tab_all <- dat2mat(dset$dat |> filter(good), "count_sam")

  map_dfr(contrasts, function(ctr) {
    conditions <- str_split(ctr, "-") |> unlist()
    meta <- dset$metadata |>
      select(sample, condition) |>
      mutate(condition = as.character(condition)) |>   # factors mess up design matrix
      filter(condition %in% conditions) |>
      distinct()
    design_mat <- model.matrix(~condition, meta)
    coef <- colnames(design_mat)[2]

    tab_all[, meta$sample] |>
      DGEList() |>
      calcNormFactors() |>
      estimateDisp(design = design_mat) |>
      glmQLFit(design = design_mat) |>
      glmQLFTest(coef = coef) |>
      topTags(n = 1e16, adjust.method = "BH", sort.by = "none") |>
      pluck("table") |>
      as_tibble(rownames = "id") |>
      mutate(contrast = ctr)
  }) |>
    drop_na() |>
    mutate(contrast = factor(contrast, levels = contrasts)) |>
    left_join(dset$info, by = "id") |>
    mutate(sig = FDR < fdr_limit & abs(logFC) > logfc_limit)
}


get_extreme_genes <- function(de, p_value_limit, logfc_limit) {
  de |>
    filter(PValue < p_value_limit & abs(logFC) > logfc_limit) |>
    pull(id) |>
    unique()
}
