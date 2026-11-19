#!/usr/bin/env Rscript
#
# Clusters the called cells of one pool on the raw counts (stage CLUSTER of
# ambient-plan), so that ambient RNA removal methods that need clusters can use them.
#
# Standard graph-based clustering with scrapper:
#   counts (Gene Expression features, called cells)
#   -> log-normalised counts (library-size factors)
#   -> top highly variable genes
#   -> PCA
#   -> shared nearest-neighbour graph
#   -> multilevel (Louvain) clustering
#
# Output: TSV with columns cell_id, cluster (one row per called cell).

suppressPackageStartupMessages({
  library(anndataR)
  library(SingleCellExperiment)
  library(data.table)
  library(scrapper)
})

# api_version >= 0.8.0 passes `--output <id>=<path>` once per declared output;
# argparser errors on unknown flags, so drop those pairs before parsing.
# (Same helper as in the other ambient-plan modules, kept here so that
# src/common/ stays an unmodified copy of the boilerplate.)
strip_output_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  keep <- rep(TRUE, length(argv))
  i <- 1
  while (i <= length(argv)) {
    if (argv[i] == "--output" && i < length(argv)) {
      keep[i] <- FALSE
      keep[i + 1] <- FALSE
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  argv[keep]
}

# arg parsing
source("src/common/cli.R")
p <- arg_parser("Clustering module (scrapper, multilevel)")
p <- add_base_args(p)
p <- add_stage_args(p, "CLUSTER")
p <- add_argument(p, "--n_hvgs", type = "integer", default = 2000, help = "Number of highly variable genes")
p <- add_argument(p, "--n_pcs", type = "integer", default = 20, help = "Number of principal components")
p <- add_argument(p, "--n_neighbors", type = "integer", default = 10, help = "Number of neighbours for the SNN graph")
p <- add_argument(p, "--resolution", type = "numeric", default = 1, help = "Multilevel (Louvain) resolution")
p <- add_argument(p, "--seed", type = "integer", default = 42, help = "Random seed (PCA and clustering)")
args <- parse_args(p, argv = strip_output_args())

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))

main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  message("  reading input files ..")
  cells_dt <- fread(args$called_cells_tsv)
  sce <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment", x_mapping = "counts")

  message("  keeping only RNA (Gene Expression) features ..")
  sce <- sce[rowData(sce)$feature_type == "Gene Expression", ]

  message(sprintf("  subsetting to %d called cells", nrow(cells_dt)))
  sce <- sce[, cells_dt$cell_id]
  counts <- counts(sce)

  message("  log-normalising ..")
  sf <- centerSizeFactors(colSums(counts))
  logcounts <- normalizeCounts(counts, sf)

  message(sprintf("  choosing %d highly variable genes ..", args$n_hvgs))
  var_stats <- modelGeneVariances(logcounts)
  hvgs <- chooseHighlyVariableGenes(var_stats$statistics$residuals, top = args$n_hvgs)

  message(sprintf("  PCA (%d PCs) ..", args$n_pcs))
  pca <- runPca(logcounts[hvgs, ], number = args$n_pcs, seed = args$seed)

  message(sprintf("  SNN graph (%d neighbours) + multilevel clustering (resolution %s) ..",
                  args$n_neighbors, args$resolution))
  graph <- buildSnnGraph(pca$components, num.neighbors = args$n_neighbors)
  clusters <- clusterGraph(graph, method = "multilevel",
                           multilevel.resolution = args$resolution, seed = args$seed)

  out <- data.table(cell_id = colnames(sce), cluster = as.character(clusters$membership))
  message(sprintf("  %d cells in %d clusters", nrow(out), length(unique(out$cluster))))
  print(table(out$cluster))

  out_path <- file.path(args$output_dir, paste0(args$name, "_clusters.tsv"))
  fwrite(out, out_path, sep = "\t")
  message(sprintf("  wrote: %s", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
