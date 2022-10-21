pivot_dat <- function(raw, suffix) {
  raw |>
    select(id = Geneid, ends_with(suffix)) |>
    pivot_longer(-id, names_to = "sample", values_to = paste0("count_", suffix)) |>
    mutate(sample = str_remove(sample, paste0(".", suffix)))
}


read_counts <- function(file, meta, cutoff = 0) {
  raw <- read_tsv(file, skip = 2, show_col_types = FALSE)

  dat <- pivot_dat(raw, "sam") |>
    full_join(pivot_dat(raw, "norm"), by = c("id", "sample"))

  # good genes: at least one condition mean > cutoff
  goods <- dat |>
    left_join(meta, by = "sample") |>
    group_by(id, condition) |>
    summarise(m = mean(count_sam, na.rm = TRUE)) |>
    ungroup() |>
    group_by(id) |>
    summarise(max_m = max(m)) |>
    mutate(good = max_m > cutoff) |>
    select(-max_m)

  dat <- dat |>
    left_join(goods, by = "id")

  info <- raw |>
    select(id = Geneid, chr = Chr, start = Start, end = End, strand = Strand, gene_symbol = Name, description = Desc) |>
    mutate(feature_id = if_else(is.na(gene_symbol), id, gene_symbol))

  list(
    metadata = meta,
    dat = dat,
    info = info
  )
}


dat2mat <- function(dat, what, id_col = "id") {
  dat |>
    rename(id = !!id_col, val = !!what) |>
    pivot_wider(id_cols = id, names_from = sample, values_from = val) |>
    column_to_rownames("id") |>
    as.matrix()
}


