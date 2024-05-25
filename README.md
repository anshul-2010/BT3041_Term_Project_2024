# A Multimodal Approach: Integrating scRNA-seq Analysis with Supervised Learning for OSF Biomarker Classification

This repository provides a comprehensive analysis of single-cell RNA-sequencing (scRNA-seq) data from Oral Submucous Fibrosis (OSF) tissue samples. The project delves into the cellular heterogeneity of OSF, aiming to:

* **Unravel Cell Diversity**: Identify and characterize distinct cell populations present within OSF tissue using unsupervised clustering techniques.
* **Cell Type Annotation**: Leverage marker gene expression and reference literature to accurately assign cell type identities to the identified clusters.
* **Differential Expression Analysis**: Uncover genes with significant expression changes between cell populations, potentially pinpointing crucial factors in OSF pathogenesis.
* **Supervised Learning for Biomarker Classification**: Develop a model capable of classifying novel genes into predefined biomarker categories within specific cell types, aiding in potential future applications for personalized medicine.

## Project Structure
The repository is organized as follows:

* `Data`: This folder houses the raw scRNA-seq count matrices for both control and experimental samples. Additionally, you'll find metadata files containing relevant information about each sample and a reference file listing established 
human cell cycle marker genes.
* `Methods `: This directory contains the codes for documenting the entire analysis workflow, including:
    * `Preprocessing and Quality Control`: Explore data cleaning steps like filtering low-quality cells and normalizing gene expression data.
    * `Cell Clustering`: Delve into the process of identifying distinct cell populations using unsupervised clustering methods like Seurat.
    * `Cell Type Annotation`: Learn how cell clusters are meticulously assigned cell type labels based on marker gene expression patterns and insights gleaned from relevant scientific literature.
    * `Differential Expression Analysis`: Discover genes exhibiting significant expression changes between cell populations, potentially offering clues into the underlying mechanisms of OSF development.
    * `Supervised Learning Model Development`: Witness the creation and training of a model designed to classify novel genes into predefined biomarker classes specific to identified cell types.
* `Images`: This folder stores all the visualizations generated throughout the analysis, including UMAP plots, violin plots, and heatmaps, offering visual representations of the data and the results.

## Running the analysis
### Environment Setup:
Install RStudio, Python and the required libraries. Tools like conda or pip can be used for dependency management.
### Data Acquisition:
Download the scRNA-seq data from a relevant database (link to drive is given)
Place the downloaded data files in the data folder, ensuring they match the file structure expected by the notebooks.
### Notebook Execution:
For the first and second parts, run the R scripts in Rstudio.
For the third part, launch a Jupyter Notebook environment. Open the notebook and execute them sequentially.

## Future Directions
* `Pseudo-Bulk RNA-Seq Analysis`: Integrate pseudo-bulk RNA-seq analysis to pinpoint genes with the most significant contributions to the OSF condition.
* `Model Optimization`: Refine the supervised learning model by optimizing hyperparameters, exploring feature engineering techniques, and incorporating additional domain knowledge to enhance its accuracy and interpretability.
* `Clinical Applications`: Investigate the potential of the developed model for personalized medicine approaches in OSF, enabling tailored treatment strategies based on individual patient profiles.
