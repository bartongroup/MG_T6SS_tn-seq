# Analysis of Ssp4 and Ssp6 effects on Pseudomonas fluorescens using Tn-seq

The R code to accompany the manuscript "Two ion-selective pore forming effectors delivered by the same Type VI secretion system have distinct consequences for targeted bacterial cells" by Mark Reglinski, Quenton W Hurst, David Williams, Marek Gierlinski, Katharine Mathers, Adam Ostrowski, Megan Bergkessel, Ulrich Zachariae, Samantha J Pitt and Sarah J Coulthurst.


## Instructions

The R code in this repository uses `renv` and `targets` for reproducibility. `renv` is used to reproduce the R environment:

```
install.packages("renv")
renv::restore()
```

Next, the [targets](https://books.ropensci.org/targets/) pipeline is invoked with the following command:

```
targets::tar_make()
```

Finally, the HTML report, including all figures, can be created with:

```
quarto::quarto_render("doc/report.qmd")
```

The HTML file will be created in folder `doc`.
