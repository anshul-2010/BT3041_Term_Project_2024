# Seurat - Guided Clustering Tutorial

## Setup the Seurat Object
For this tutorial, we will be analyzing a dataset of Oral submucous fibrosis (OSF), which was made available to public on May 17, 2024. It is freely available from 10X Genomics. There are 9,124 single cells that were sequenced.

We start by reading in the data. The `Read10X` function reads in the output of the **cellranger** pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e., gene: row) that are detected in each cell (column). Note that more recent versions of cellranger now also has the functionality to output h5 file format.

WE next use the count matrix to create a Seurat object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset. For example, in Seurat v5, the count matrix is stored in **pbmc[["RNA"]]$counts**.

On printing the count matrix, the . values in the matrix represent 0s (no molecules detected). Since most values in an scRNA-seq matrix are 0, Seurat uses a sparse-matrix representation whenever possible. This results in significant memory and speed savings for Drop-seq/inDrop/10x data.

## Standard pre-processing workflow
The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. These represent the selection and filtration of cells based on QC metrics, data normalization and scaling, and the detection of highly variable features.

### QC and selecting cells for further analysis
Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. A few QC metrics commonly used by the community include

* The number of unique genes detected in each cell:
    * Low-quality cells or empty droplets will often have very few genes.
    * Cell doublets or multiplets may exhibit an aberrantly high gene count.
* Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes).
* The percentage of reads that map to the mitochondrial genome:
    *Low-quality / dying cells often exhibit extensive mitochondrial contamination.
    * We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features.
    * We use the set of all genes starting with `MT-` as a set of mitochondrial genes.
* The number of unique genes and total molecules are automatically calculated during `CreateSeuratObject()`.
    * You can find them stored in the object meta data.

<img src="https://github.com/anshul-2010/Computational-Systems-Biology/blob/main/images/display/Actor_Critic.png" alt="Actor Critic model" width="350"/>



<img src="https://github.com/anshul-2010/Computational-Systems-Biology/blob/main/images/display/RL-dFBA.png" alt="RL-dFBA working" width="500"/>