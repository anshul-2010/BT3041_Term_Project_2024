# Differential Expression Analysis of Integrated scRNA-seq Data
This section details the integration and differential expression (DE) analysis pipeline employed to identify genes with expression changes associated with the experimental condition of interest in OSF.

## Data Integration and Quality Control
* **Control and Experimental Datasets**: scRNA-seq data from both control and experimental groups are integrated for robust DE analysis.
* **Quality Control (QC)**: A comprehensive QC step ensures data integrity. Metrics like UMI counts, detected genes, mitochondrial expression, and cell complexity are assessed to filter low-quality cells.
* **Visualization**: Pre- and post-QC distributions of these metrics are visualized to demonstrate the effectiveness of filtering.

## Preprocessing
* **SCTransform**: Seurat's `SCTransform` function performs normalization, variance stabilization, and covariate regression to improve expression estimates.
* **Cell Cycle Correction**: Cell cycle variation is accounted for using Seurat's `CellCycleScoring` function and VariableFeatures selection.
* **Dimensionality Reduction**: Principal Component Analysis (PCA) explores major sources of variation within the preprocessed data.

## Integration Strategy
* **Canonical Correlation Analysis (CCA)**: CCA identifies shared variation between datasets for initial alignment.
* **Mutual Nearest Neighbors (MNNs)**: MNNs across datasets serve as anchors for linking and integrating the data.
* **Visualization**: UMAPs visualize the data before and after integration.

## Differential Expression Analysis Workflow
* **Graph-based Clustering**: Clustering identifies shared condition-specific cell populations based on integrated data.
* *Validation*: Distribution of cells across samples and clusters is examined to assess potential batch effects.
* **Cell Cycle and Mitochondrial Effects**: Cell cycle phase and mitochondrial expression are visualized across clusters to ensure minimal confounding effects.
* **Metrics Exploration**: UMI counts, genes per cell, and separation based on PCs are evaluated.

## Marker Gene Identification
* **Differential Expression (DE) Analysis**: DE analysis identifies genes significantly overexpressed within specific clusters compared to others.
* **Marker Gene Visualization**: Violin plots visualize expression patterns of identified marker genes (Figures 41-45).

## Key Findings
* Clusters 1, 2, and 7 exhibit consistent marker gene expression, suggesting a potential lack of direct association with OSF pathogenesis.
* Clusters 0 and 5 display significant changes in gene expression, highlighting genes potentially crucial for OSF development or progression.
* Genes like AGR2 and FAM3B exhibit altered cluster membership, potentially indicating involvement in OSF.

## Limitations and Future Work
* Current DE analysis methods might overestimate significance due to treating individual cells as independent samples.
* Future work could explore incorporating pseudo-bulk DE analysis workflows for more accurate conclusions.

This README provides a high-level overview of the differential expression analysis section.  Remember to replace the bracketed figure references with the actual figure names or links within your repository.