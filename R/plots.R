okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#AAAAAA", "#000000")


plot_count_distribution <- function(dset, cutoff = 1) {
  dset$dat |>
    ggplot(aes(x = count_sam, y = ..density..)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_histogram(bins = 60) +
    geom_vline(xintercept = cutoff, colour = "red", linetype = "dashed") +
    facet_wrap(~sample) +
    scale_x_log10() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Count", y = "Density")
}

plot_genes <- function(dset, gids, what = "count_norm", point_size = 3) {
  dset$dat |>
    filter(id %in% gids) |>
    left_join(dset$metadata, by = "sample") |>
    left_join(dset$info, by = "id") |>
    mutate(
      condition = factor(condition, levels = levels(dset$metadata$condition)),
      gene_symbol = if_else(is.na(gene_symbol), id, gene_symbol)
    ) |>
  ggplot(aes(x = condition, y = get(what), fill = replicate)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_point(shape = 21, colour = "grey50", size = point_size, position = position_dodge(width = 0.2)) +
    scale_fill_manual(values = okabe_ito_palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    facet_wrap(~ gene_symbol, scales = "free_y") +
    labs(x = NULL, y = what)
}


plot_replicate_comparison <- function(dset, what = "count_norm") {
  tab <- dat2mat(dset$dat, what)
  xydat <- dset$metadata |>
    group_split(condition) |>
    map_dfr(function(w) {
      expand_grid(x = w$sample, y = w$sample) |>
        filter(x < y) |>
        unite(pair, c(x, y), sep = ":", remove = FALSE)
    }) |>
    rowwise() |>
    group_split() |>
    map_dfr(function(r) {
      tibble(
        pair = r$pair,
        id = rownames(tab),
        x = tab[, r$x],
        y = tab[, r$y]
      )
    })

  xydat |>
    ggplot(aes(x = x, y = y)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_point(size = 0.1) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    facet_wrap(~pair) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = NULL, y = NULL)
}


plot_volma <- function(res, p, fdr, fc, group, point_size, point_alpha) {
  r <- res |>
    mutate(
      p = get(p),
      group = get(group)
    ) |>
    select(x, y, sig, group)
  rsel <- r |> filter(sig)
  rm(res)  # Minimise environment for serialisation
  ggplot(r, aes(x = x, y = y)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(size = point_size, alpha = point_alpha, colour = "grey50") +
    geom_point(data = rsel, colour = "black") +
    facet_grid(. ~ group)
}

plot_ma <- function(res, a = "logCPM", fc = "logFC", p = "PValue", fdr = "FDR", group = "contrast",
                    point_size = 0.5, point_alpha = 0.5) {
  res |>
    mutate(
      x = get(a),
      y = get(fc),
    ) |>
    plot_volma(p, fdr, fc, group, point_size, point_alpha) +
    geom_hline(yintercept = 0, size = 0.1, alpha = 0.5) +
    labs(x = expression(log[10]~Intensity), y = expression(log[2]~FC))
}

plot_volcano <- function(res, fc = "logFC", p = "PValue", fdr = "FDR", group = "contrast",
                         point_size = 0.5, point_alpha = 0.5) {
  res |>
    mutate(
      x = get(fc),
      y = -log10(get(p)),
    ) |>
    plot_volma(p, fdr, fc, group, point_size, point_alpha) +
    geom_vline(xintercept = 0, size = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
}

