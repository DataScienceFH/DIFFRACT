![DIFFRACT](DIFFRACT.jpg)

# DIFFRACT: minimal working example

DIFFRACT quantifies condition-dependent changes in molecular network geometry
using local fractal dimension (LFD). This repository provides a minimal working
R implementation of its core numerical workflow and a compact example dataset.
It demonstrates edge weighting, node-wise LFD calculation, structural/dynamic
decomposition, comparison of two conditions and extraction of induced components.

The implementation accompanies the preprint
[*Local fractal dimension as an indicator of functional reorganization in
molecular interaction networks*](https://doi.org/10.64898/2026.10.27.762149).

## Repository contents

| File | Description |
|---|---|
| `R/DIFFRACT.R` | Complete, executable analysis script |
| `example_data/STRINGfull_v12_sample.graphml` | Connected, undirected network with 1,000 nodes and 999 edges |
| `example_data/proteinabundance_WT_KO.csv` | Matching abundance matrix with three WT and three KO samples |

The script and the two example-data files are sufficient to run the analysis.
Results are generated in an `output/` directory beside the script.

## Installation and usage

The example was tested with **R 4.5.2**, **igraph 2.2.3**, **infotheo 1.2.0.1**
and **energy 1.7-12**. Install the dependencies once:

```r
utils::install.packages(c("igraph", "infotheo", "energy"))
```

Download the repository and retain the file layout shown above. In RStudio,
open `DIFFRACT.R` and click **Source**. Alternatively, from R with the
repository folder as the working directory:

```r
base::source("DIFFRACT.R")
```

From a terminal in the repository folder:

```text
Rscript DIFFRACT.R
```

The script resolves input and output paths relative to its own location and
processes all four weighting strategies automatically. After package
installation, the example runs locally without downloads or additional data.

## Example data

### Network

The network is a connected induced subset of the human **STRING v12.0**
protein-association network (von Mering et al., 2003; Szklarczyk et al., 2025;
see [References](#references)), restricted to proteins represented in the
abundance dataset. Associations were retained at a combined confidence score
of at least **900** on the 0-1,000 scale.

The GraphML vertex attribute `name` contains neutral IDs **N0001-N1000**.
These match the abundance table exactly.

### Abundance matrix

The example is derived from proteomics with three wild-type and three
knockout samples, retaining proteins with complete measurements and
positive variance across the six samples.

| Column | Contents |
|---|---|
| `node` | Neutral ID matching the GraphML vertex attribute `name` |
| `WT1`, `WT2`, `WT3` | Transformed wild-type abundances |
| `KO1`, `KO2`, `KO3` | Transformed knockout abundances |

## Analysis workflow

### 1. Condition-specific edge weighting

WT and KO profiles are analyzed separately on the same fixed scaffold.

| Strategy | Association or control metric |
|---|---|
| **MI** | Mutual information after equal-frequency discretization into two bins |
| **Pearson** | Absolute Pearson correlation |
| **dCor** | Distance correlation calculated with `energy::dcor` |
| **VarMiss** | Constant divided by `((var(x) + var(y))/2 + epsilon)^beta`, with `epsilon = 0.01` and `beta = 4` |

### 2. Node-wise LFD and decomposition

For each node, shortest-path distances determine the cumulative number of
other nodes reached at each distinct positive radius. The slope of
`log(cumulative node count) ~ log(radius)` is the LFD. All positive finite
radii are used, with at least two radii required for fitting. The focal node
is excluded from the count; regression R-squared is saved with each estimate.

Structural LFD uses unweighted shortest paths and is calculated once.
Weighted LFD is calculated independently for WT and KO under each strategy:

```text
LFD_dynamic = LFD_weighted - LFD_structural
delta_LFD = LFD_dynamic(KO) - LFD_dynamic(WT)
```

The dynamic contribution is stored as `lfd_modulated` in the output tables.
Positive `delta_lfd` indicates higher dynamic LFD in KO than in WT.

### 3. Differential LFD and induced components

Node-wise differences are standardized across the network for each strategy.
Nodes with **p < 0.05** define an induced subgraph on the fixed scaffold.
Both directions of LFD change enter the same graph. Connected components
with **at least five nodes** are exported individually. The script checks
that each strategy produces at least one qualifying component; replacement
data that fail this demonstration requirement trigger an error after the
analysis results have been written.

## Generated outputs

Paths below are relative to `output/`. `<method>` denotes `MI`, `Pearson`,
`dCor` or `VarMiss`.

| Output | Contents |
|---|---|
| `LFD_structural.csv` | Unweighted LFD (`nfd`) and regression R-squared (`r2`) for each node |
| `parameters.csv` | Recorded weighting, distance, contrast and extraction parameters |
| `sessionInfo.txt` | R version and package versions for the run |
| `component_summary.csv` | Selected-node, selected-edge and retained-component counts |
| `<method>/LFD_WT.csv`, `<method>/LFD_KO.csv` | Weighted `nfd`, `r2`, `lfd_structural` and `lfd_modulated` |
| `<method>/delta_LFD_KO_vs_WT.csv` | WT and KO dynamic LFD, `delta_lfd`, `z_score` and `p_value` for all nodes |
| `<method>/weighted_WT.graphml`, `<method>/weighted_KO.graphml` | Condition-specific networks with edge distances in `weight` |
| `<method>/components/significant.graphml` | Induced graph of all selected nodes, including isolated nodes |
| `<method>/components/membership.csv` | Component IDs, sizes, retention flags and exported component numbers |
| `<method>/components/GIANT.graphml` | Largest selected-node component, exported when the selected graph has edges |
| `<method>/components/Component_<n>.graphml` | Each retained component with at least five nodes |
| `<method>/components/Component_<n>.csv` | Differential LFD table for the corresponding retained component |

Component graphs include node attributes `p_value` and `deltaLFD_modulated`.
In `membership.csv`, `exported_component` links to the `Component_<n>` files
and is `NA` for discarded components. All node-level outputs use neutral IDs.

### Expected component extraction

With the provided inputs and tested environment:

| Strategy | Selected nodes | Edges between selected nodes | Retained components | Retained component size |
|---|---:|---:|---:|---:|
| MI | 42 | 6 | 1 | 5 |
| Pearson | 55 | 9 | 1 | 6 |
| dCor | 55 | 11 | 1 | 6 |
| VarMiss | 23 | 14 | 1 | 13 |

## Citation

Please cite the accompanying DIFFRACT preprint when using this implementation:

Hause, F. *et al.* Local fractal dimension as an indicator of functional
reorganization in molecular interaction networks. *bioRxiv* (2026).
[https://doi.org/10.64898/2026.10.27.762149](https://doi.org/10.64898/2026.10.27.762149).

## References

The following STRING references are cited in the accompanying manuscript:

1. von Mering, C. *et al.* STRING: a database of predicted functional
   associations between proteins. *Nucleic Acids Res.* **31**, 258-261 (2003).
2. Szklarczyk, D. *et al.* The STRING database in 2025: protein networks with
   directionality of regulation. *Nucleic Acids Res.* **53**, D730-D737 (2025).

## Contact

**Frank Hause**  
Email: [**frank.hause@medizin.uni-halle.de**](mailto:frank.hause@medizin.uni-halle.de)
