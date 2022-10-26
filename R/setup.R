SAMPLES <- c("YL37_1", "YL37_2", "YL37_3", "GM103_1", "GM103_2", "GM103_3",
             "YL57_1", "YL57_2", "YL57_3")

COUNT_FILE <- "data/220707_counts.txt"
EGGNOG_FILE <- "data/eggnog/MM_clhlot17.emapper.annotations.tsv"

CONTRASTS <- c("YL57-YL37", "GM103-YL37", "GM103-YL57")

FDR_LIMIT <- 0.05
LOGFC_LIMIT <- 0

CV_LIMITS <- tibble::tribble(
  ~LIMIT,
  100,
  1,
  0.5
)

make_metadata <- function() {
  tibble(
    sample = SAMPLES
  ) |>
    separate(sample, c("condition", "replicate"), sep = "_", remove = FALSE) |>
    mutate(condition = as_factor(condition) |> fct_relevel("YL37", "YL57"))
}
