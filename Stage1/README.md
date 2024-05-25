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
    * Low-quality / dying cells often exhibit extensive mitochondrial contamination.
    * We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features.
    * We use the set of all genes starting with `MT-` as a set of mitochondrial genes.
* The number of unique genes and total molecules are automatically calculated during `CreateSeuratObject()`.
    * You can find them stored in the object meta data.

In the example below, we visualize QC metrics, and use these to filter cells.

* We filter cells that have unique feature counts over 2,500 or less than 200
* We filter cells that have >5% mitochondrial counts

## Normalizing the data
After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. In Seurat v5, Normalized values are stored in **pbmc[["RNA"]]$data**.

For clarity, in this previous line of code (and in future commands), we provide the default values for certain parameters in the function call.

While this method of normalization is standard and widely used in scRNA-seq analysis, global-scaling relies on an assumption that each cell originally contains the same number of RNA molecules. The use of SCTransform replaces the need to run NormalizeData, FindVariableFeatures, or ScaleData (described below.)

## Identification of highly variable features 
We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). Our procedure directly models the mean-variance relationship inherent in single-cell data, and is implemented in the FindVariableFeatures() function. By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.

## Scaling the data 
Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:

* Shifts the expression of each gene, so that the mean expression across cells is 0
* Scales the expression of each gene, so that the variance across cells is 1
    * This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
    * The results of this are stored in pbmc[["RNA"]]$scale.data
* By default, only variable features are scaled.
* You can specify the features argument to scale additional features.

In Seurat, we also use the ScaleData() function to remove unwanted sources of variation from a single-cell dataset. For example, we could `regress out` heterogeneity associated with (for example) cell cycle stage, or mitochondrial contamination. However, particularly for advanced users who would like to use this functionality, we strongly recommend the use of our new normalization workflow, SCTransform().

## Perform linear dimensional reduction 
Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset (if you do want to use a custom subset of features, make sure you pass these to ScaleData first).

For the first principal components, Seurat outputs a list of genes with the most positive and negative loadings, representing modules of genes that exhibit either correlation (or anti-correlation) across single-cells in the dataset. Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap().

In particular DimHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and features are ordered according to their PCA scores. Setting cells to a number plots the `extreme` cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated feature sets.

## Determine the dimensionality of dataset 
To overcome the extensive technical noise in any single feature for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a `metafeature` that combines information across a correlated feature set. The top principal components therefore represent a robust compression of the dataset. An alternative heuristic method generates an `Elbow plot`: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot() function). In this example, we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.

## Run non-linear dimensionality reduction (UMAP)
Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. The goal of these algorithms is to learn underlying structure in the dataset, in order to place similar cells together in low-dimensional space. Therefore, cells that are grouped together within graph-based clusters determined above should co-localize on these dimension reduction plots.

While we and others have routinely found 2D visualization techniques like tSNE and UMAP to be valuable tools for exploring datasets, all visualization techniques have limitations, and cannot fully represent the complexity of the underlying data. In particular, these methods aim to preserve local distances in the dataset (i.e. ensuring that cells with very similar gene expression profiles co-localize), but often do not preserve more global relationships. We encourage users to leverage techniques like UMAP for visualization, but to avoid drawing biological conclusions solely on the basis of visualization techniques.

## Finding differentially expresed features (cluster biomarkers)
Seurat can help you find markers that define clusters via differential expression (DE). By default, it identifies positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. FindAllMarkers() automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

In Seurat v5, we use the presto package (as described here and available for installation here), to dramatically improve the speed of DE analysis, particularly for large datasets. For users who are not using presto, you can examine the documentation for this function (?FindMarkers) to explore the min.pct and logfc.threshold parameters, which can be increased in order to increase the speed of DE testing.

We include several tools for visualizing marker expression. VlnPlot() (shows expression probability distributions across clusters), and FeaturePlot() (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations. DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster. Fortunately in the case of this dataset, we can use canonical markers to easily match the unbiased clustering to known cell types.
