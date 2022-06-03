---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# BioLiP clusters

On this page, I'll apply the cluster labels to the BioLiP data, and investigate
the distribution of clusters.


## Setup

Add user's current working directory to search path:
```{code-cell}
import sys
sys.path = [r"C:\sw\plb"] + sys.path
```

Also change current working directory:
```{code-cell}
import os
os.chdir(r"C:\sw\plb")
```

Imports:
```{code-cell}
from collections import Counter
from pathlib import Path

import bokeh.io
import myst_nb

from src.config import get_config
from src.data_structure_utils import read_tsv
from src.stats_utils import (
    get_stats_from_counter,
    print_counter_stats_df,
)
from src.plotting_functions import (
    plot_cluster_counts_bdb_biolip_combined,
    plot_cluster_counts_single_dataset,
)

```

Set up document options:
```{code-cell}
# This import automatically sets seeds
import src.set_seed

# Nice pandas tables
import itables
itables.init_notebook_mode(all_interactive=True)
itables.options.maxColumns = 0  # unlimited columns
itables.options.maxBytes = 32000000
itables.options.classes = ["display", "nowrap"]
itables.options.order = [] # disable auto-sorting

# Set up bokeh interactive plots
bokeh.io.output_notebook()
```

Get filepaths of data directories:
```{code-cell}
config = get_config()
dir_data = Path(config["data_directory"])
dir_features = Path(config["features_directory"])
```

Load dataframe with cluster assignments:
```{code-cell}
df_biolip_ann_with_clusters = read_tsv(dir_features / "df_biolip_ann_with_clusters.tsv")
```

## Distribution of protein clusters

Similarly to BindingDB, we can generate statistics for the protein cluster
counts in BioLiP:

```{code-cell}
biolip_cluster_counts_p: Counter[str] = Counter(
    df_biolip_ann_with_clusters["protein_cluster"]
)
biolip_cluster_stats_p = get_stats_from_counter(
    counts=biolip_cluster_counts_p,
    ntop=30,
    nbins=40,
    name_of_thing_being_counted="protein_cluster_counts",
    xaxis_name="number of occurrence of protein cluster",
    do_show=True,
)
```

Most commonly occuring protein clusters:
```{code-cell}
:tags: ["remove-input"]
print_counter_stats_df(biolip_cluster_stats_p)
```


## Distribution of ligand clusters

Let's do the same, but for the ligand clusters:

```{code-cell}
biolip_cluster_counts_l: Counter[str] = Counter(
    df_biolip_ann_with_clusters["ligand_cluster"]
)
biolip_cluster_stats_l = get_stats_from_counter(
    counts=biolip_cluster_counts_l,
    ntop=30,
    nbins=40,
    name_of_thing_being_counted="ligand_cluster_counts",
    xaxis_name="number of occurrences of ligand cluster",
    do_show=True,
)
```

Most commonly occuring ligand clusters:
```{code-cell}
:tags: ["remove-input"]
print_counter_stats_df(biolip_cluster_stats_l)
```

## Distribution of (p, l) cluster tuples

Finally, we can show the statistics and plots for the $(p,l)$ cluster tuples
present in the BioLiP data like we did for the BindingDB data:

```{code-cell}
biolip_cluster_counts_pl: Counter[tuple[str, str]] = Counter(
    zip(
        df_biolip_ann_with_clusters["protein_cluster"],
        df_biolip_ann_with_clusters["ligand_cluster"],
    )
)
biolip_cluster_stats_pl = get_stats_from_counter(
    counts=biolip_cluster_counts_pl,
    ntop=30,
    nbins=40,
    name_of_thing_being_counted="pl_cluster_tuple",
    xaxis_name="number of occurrences of (p, l) cluster tuple",
    do_show=True,
)
```

The most commonly occuring $(p,l)$ cluster pairs are:
```{code-cell}
print_counter_stats_df(biolip_cluster_stats_pl)
```

Heatmap (after reindexing $p$ and $l$ values):

```{code-cell}
plot_cluster_counts_single_dataset(
    pl_counts=biolip_cluster_counts_pl, backend="bokeh", permute_rows_and_columns=False
)
```

The same heatmap, but with permuted rows and columns (to put most of the points
near the top left):

```{code-cell}
plot_cluster_counts_single_dataset(
    pl_counts=biolip_cluster_counts_pl, backend="bokeh", permute_rows_and_columns=True
)
```
