BOX# CellSP
__Note: Repository is work in progress__
<br>
<br>
CellSP is a python package for the analysis of subcellular spatial transcriptomic data. CellSP works with datasets generated at single-modulecule resolution from technologies like Xenium, CosMx, MERSCOPE or other ISH-like data. Using existing tools [InSTAnT](https://github.com/bhavaygg/InSTAnT) and [SPRAWL](https://github.com/salzman-lab/SPRAWL/), CellSP identifies statistically signficant subcellular patterns of gene transcripts and uses a biclustering algorithm to aggregate these patterns over hundereds of cells to produce "gene-cell modules". These modules represent the consistent detection of the same subcellular pattern by a set of genes in the same cells and offer a summarized and biologically interpretable desciption of subcellular patterns. CellSP provides specialized techniques for visualizing such modules and their defining spatial patterns. Additionally, CellSP utilize Gene Ontology (GO) enrichments tests to offer functionsal insights into the genes comprising the module as CellSPll as the cells comprising the module.

![CellSP_overview](https://github.com/bhavaygg/CellSP/blob/main/figures/Overview.png)

***

## How to install CellSP

CellSP recommend using our environment.yml file to create a new conda environment to avoid issues with package incompatibility.

```
conda env create -f environment.yml
```
This will create a new conda environment with the name `CellSP` and has all dependencies installed. 

Alternatively, the package can be installed using pip.

```
pip install cellSP
```

***
## How to use CellSP

CellSP expects data (both single cell and spatial transcriptomic) to be in AnnData format and can be loaded using 

```
adata_sc, adata_st = cellSP.ds.load_data(sc_adata= 'files/adata_sc.h5ad', st_adata = "files/adata_st.h5ad")
```

**Note - Single cell data on the same tissue is required for characterization of the module cells.**

To load raw csv data, refer to [file]() for instructions.

CellSP preprocess the input single cell data by performing denoising using [MAGIC](https://github.com/KrishnaswamyLab/MAGIC) and impute the expression of genes not in the ST panel using [Tangram](https://github.com/broadinstitute/Tangram/). 

```
adata_sc = cellSP.pp.impute(adata_sc, t="auto")
adata_st = cellSP.pp.run_tangram(adata_sc, adata_st, device='cpu')
```

After Tangram imputation, the single cell and spatial Anndata objects are combined into one. This completes the preprocessing required for using CellSP. This can be skipped if cellular characterization is not required.

There are three main steps involved in running CellSP - 
1.  Subcellular Pattern Discovery
2.  Module Discovery
3.  Module Characterization

### Subcellular Pattern Discovery

CellSP uses InSTAnT and SPRAWL for identifying statistically significant subcellular patterns. InSTAnT tests if transcripts of a gene pair tend to be proximal to each other more often than expected by chance, while SPRAWL identifies four types of subcellular patterns – peripheral, radial, punctate and central – describing the distribution of a gene’s transcripts within the cell. 

To run InSTAnT, CellSP has two primary parameters - 
    - `distance_threshold`: The distance (in microns) at which to consider 2 genes proximal.
    - `alpha_cpb`: p-value signifiance threshold below which a gene-pair is considered colocalized for the CPB test. Default = 1e-3

```
adata_st = cellSP.ch.run_instant(adata_st = adata_st, distance_threshold=2, alpha_cpb=1e-5)
```

To run SPRAWL, CellSP uses the default parameters from the original implementation.

```
adata_st = cellSP.ch.run_sprawl(adata_st)
```

### Module Discovery

CellSP use a biclustering tool, LAS, to analyze each of the patterns and identiy "gene-cell modules". Each module represents a set of genes or gene pairs that exhibit the same type of sub-cellular pattern in the same set of cells, with statistical significance estimated by a Bonferroni-based score.

CellSP has 2 functions for module discovery, one for SPRAWL and one for InSTAnT. Both the functions share the same parameters but the InSTAnt function has two additional parameter 
    - `alpha`: p-value signifiance threshold below which a gene-pair is considered for biclustering. Default = 1e-3
    - `topk`: Select only the K most significant gene pairs that have p-value < `alpha`. Default = None

These parameters is used the restrict the number of gene-pairs over which biclustering is performed in order to reduce the computational complexity. 

The other parameters used are - 
    - `num_biclusters`: Number of modules to find. Default = 10.
    - `randomized_searches`: Number of randomized searches to perform in LAS. Default = 50000.

```
adata_st = cellSP.ch.bicluster_instant(adata_st, distance_threshold=2, threads=128, alpha=1e-5, num_biclusters = 50, randomized_searches = 50000)
adata_st = cellSP.ch.bicluster_sprawl(adata_st, threads=128, num_biclusters = 50, randomized_searches = 50000)
```

### Module Characterization

To aid biological interpretation, CellSP reports shared properties of the genes and cells of each discovered module. Genes are characterized using Gene Ontology (GO) enrichment tests, while cells are characterized by their cell type composition if such information is available. To provide a more precise characterization of a module’s cells, CellSP trains a machine learning classifier to discriminate those cells from all other cells, using the expression levels of all genes other than the module genes. Genes that are highly predictive in this task are then subjected to GO enrichment tests, furnishing hypotheses about biological processes and pathways that are active specifically in the module cells.

To characterize the module genes -

```
adata_st = cellSP.geo.geo_analysis(adata_st, setting="module")
```


To characterize the module cells, we first train a random forest classifier to find genes that are predictive of module presence and then perform enrichment tests -  
```
adata_st = cellSP.md.model_modules(adata_st, do_shap=True, subsample = True)
adata_st = cellSP.geo.geo_analysis(adata_st, setting="cell")
```

### Visualization

To help visualize modules defined by the five types of subcellular spatial patterns (four types identified by SPRAWL and colocalization patterns identified by InSTAnT), we developed three complementary plotting techniques.

![CellSP_visualizations](https://github.com/bhavaygg/CellSP/blob/main/figures/CellSP_visualizations.png)


SPRAWL identifies spatial localization patterns (peripheral, central, radial, or punctate) for each gene in individual cells. CellSP aggregates these patterns across cells by representing each cell within a standardized unit circle, enabling comparative analyses.

For “central” and “peripheral” patterns, gene densities are averaged across concentric rings of the circle, revealing expected trends: higher densities in the innermost rings for central patterns and in the outermost rings for peripheral patterns.

For “punctate” and “radial” patterns, the circle is divided into sectors, and densities are aligned to highlight directional concentration. Module genes are expected to cluster in specific sectors, while non-module genes display uniform distributions.

InSTAnT identifies colocalized gene pairs based on their spatial proximity. A proximity enrichment score compares the colocalization of gene pairs in module cells versus non-module cells. Visualization includes heatmaps showing enrichment scores for module and control genes, highlighting patterns of colocalization and contrasting them with non-module gene behavior.
***

### How to cite CellSP

```
@article{aggarwal2025cellsp,
  title={CellSP: Module discovery and visualization for subcellular spatial transcriptomics data},
  author={Aggarwal, Bhavay and Sinha, Saurabh},
  journal={bioRxiv},
  pages={2025--01},
  year={2025},
  publisher={Cold Spring Harbor Laboratory}
}
```
