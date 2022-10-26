pivot_dat <- function(raw, suffix) {
  raw |>
    select(id = Geneid, ends_with(suffix)) |>
    pivot_longer(-id, names_to = "sample", values_to = paste0("count_", suffix)) |>
    mutate(sample = str_remove(sample, paste0(".", suffix)))
}

remove_zero_genes <- function(dat) {
  dat |>
    group_by(id) |>
    mutate(mx = max(count_sam)) |>
    filter(mx > 0) |>
    ungroup() |>
    select(-mx)
}


read_counts <- function(file, meta) {
  raw <- read_tsv(file, skip = 2, show_col_types = FALSE)

  dat <- pivot_dat(raw, "sam") |>
    full_join(pivot_dat(raw, "norm"), by = c("id", "sample")) |>
    remove_zero_genes()

  sel <- dat$id |>
    unique()

  ms <- dat |>
    left_join(meta, by = "sample") |>
    group_by(id, condition) |>
    summarise(m = mean(count_norm), s = sd(count_norm)) |>
    ungroup() |>
    mutate(cv = s / m)

  info <- raw |>
    select(id = Geneid, chr = Chr, start = Start, end = End, strand = Strand, gene_symbol = Name, description = Desc) |>
    mutate(feature_id = if_else(is.na(gene_symbol), id, gene_symbol))

  list(
    metadata = meta,
    dat = dat,
    sel = sel,
    info = info,
    ms = ms
  )
}


dat2mat <- function(dat, what, id_col = "id") {
  dat |>
    rename(id = !!id_col, val = !!what) |>
    pivot_wider(id_cols = id, names_from = sample, values_from = val) |>
    column_to_rownames("id") |>
    as.matrix()
}

filter_cv <- function(dset, cv_limit) {
  bads <- dset$ms |>
    filter(is.na(cv) | cv > cv_limit) |>
    pull(id) |>
    unique()
  goods <- setdiff(unique(dset$dat$id), bads)
  dset$sel <- goods
  dset
}
