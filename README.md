# ambient-clustering

Clustering of called cells on the raw counts, for the ambient RNA benchmark ([ambient-plan](https://github.com/omni-scrna/ambient-plan)).

An [omnibenchmark](https://omnibenchmark.org) module for the stage `CLUSTER`, which runs after cell calling (`CELLCALL`), once per pool. The clusters are meant as input for ambient RNA removal methods that need them.

## Entrypoint `default` (`run.R`)

Standard graph-based clustering with [scrapper](https://github.com/libscran/scrapper): Gene Expression features of the called cells → log-normalised counts → highly variable genes → PCA → shared nearest-neighbour graph → multilevel (Louvain) clustering.

| flag | meaning |
|---|---|
| `--rawdata_h5ad` | raw droplets of one pool (h5ad, X = counts) |
| `--called_cells_tsv` | called barcodes (column `cell_id`) |
| `--output_dir`, `--name` | set by omnibenchmark |
| `--n_hvgs` (2000), `--n_pcs` (20), `--n_neighbors` (10), `--resolution` (1), `--seed` (42) | method parameters |

Output: `<name>_clusters.tsv` with columns `cell_id`, `cluster`.

Environment: the plan's `r_conda` env (`envs/r.yml`: scrapper, anndataR, argparser, data.table).

`src/common/` is copied from [boilerplate](https://github.com/omni-scrna/boilerplate) with `scripts/pull.py`; don't edit it by hand. The stage schema `src/common/schema/CLUSTER.json` is written here because ambient-plan has no `schema/` folder yet.
