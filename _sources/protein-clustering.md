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

# Protein clustering process

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

import myst_nb

from src.config import get_config
from src.data_split.protein_clustering import read_protein_cluster_labels
from src.data_structure_utils import read_tsv
from src.stats_utils import (
    get_stats_from_counter,
    print_counter_stats_df,
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
```

Get filepaths of data directories:
```{code-cell}
config = get_config()
dir_data = Path(config["data_directory"])
dir_features = Path(config["features_directory"])
```

## Overview

The general process I used for clustering the protein sequences is:

- Get a combined list of all unique protein sequences present in the filtered
  BindingDB and BioLiP dataframes
- Cluster the sequences using [mmseqs cluster](https://github.com/soedinglab/MMseqs2)
  with a sequence identity threshold of >10% (otherwise, I used default options)


## Cluster statistics

Read the protein cluster assignments file:
```{code-cell}
df_protein_clusters = read_protein_cluster_labels(
    proteins_fasta=dir_features / "all_proteins_single.fasta",
    clusters_tsv=dir_features / "all_proteins_single_mmseqs_db_seqid_0.1/clusters.tsv",
)
```

Print a sample of the clustering outputs:
```{code-cell}
itables.show(df_protein_clusters.head(30))
```

Generate some statistics from the clustering outputs:
```{code-cell}
protein_cluster_stats = get_stats_from_counter(
    counts=Counter(df_protein_clusters["cluster"]),
    ntop=30,
    nbins=40,
    name_of_thing_being_counted="cluster_label",
    xaxis_name="number of unique aa sequences in the cluster",
    do_show=True,
)
```

The most commonly occurring clusters are:
```{code-cell}
print_counter_stats_df(protein_cluster_stats)
```
