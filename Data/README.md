The control and experimental data for OSF condition is given in the drive link
https://drive.google.com/drive/folders/1pJ7Z5TTrfHf_qF-I3kOBx31L-j0BbcJb?usp=sharing

# Oral Submucous Fibrosis (ORF) scRNA-seq data
This repository contains single-cell RNA sequencing (scRNA-seq) data from a study investigating the cellular mechanisms underlying Oral Submucous Fibrosis (OSF). The goal was to identify novel therapeutic targets and explore the potential of Dental Pulp Stem Cells (DPSCs) as a treatment approach.

## Data Description
The data directory contains the following files:
* `expression_matrix.mtx`: This file stores the gene expression counts for each cell in the scRNA-seq experiment.
* `metadata.csv`: This CSV file provides metadata associated with each cell, including:
    * `Sample ID`: Identifier for the OSF patient or control sample.
    * `Cell Type (optional)`: Predicted cell type based on reference markers (may require further analysis)
    * `Other Features`: Additional experimental features relevant to the study design.

## Overall Design
The study employed scRNA-seq analysis to compare oral mucosal tissues from OSF patients with control samples. This analysis aimed to:

* Identify distinct cell populations within the OSF tissue.
* Investigate the functional roles of these cell populations in OSF pathogenesis.
* Characterize the communication between epithelial and immune cells in the context of OSF.

## Usage Notes
This data can be used for various downstream analyses, including:

* Differential expression analysis to identify genes differentially expressed between OSF and control samples.
* Cell type clustering to define distinct cell populations within the scRNA-seq data.
* Trajectory inference to explore potential cell fate transitions in OSF.
* Cell-cell interaction analysis to understand communication patterns between different cell types.