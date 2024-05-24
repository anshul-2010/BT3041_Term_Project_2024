# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Show the QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

##                      orig.ident nCount_RNA nFeature_RNA percent.mt
## AAACCCAAGGTAAGGA-1     pbmc3k      13920         3252  26.918103
## AAACCCAAGGTCCTGC-1     pbmc3k      12524         2984  11.609709
## AAACCCACATGGGTCC-1     pbmc3k       9196         2393   5.491518
## AAACCCAGTCACCGCA-1     pbmc3k      31777         5208   7.077446
## AAACCCAGTCTGTAGT-1     pbmc3k       9006         2284   4.696869

# Visualize QC metrics as a violin plot
violin_plot <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\violinplot.jpeg", width=600, height=350)
violin_plot
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e., columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
feature_scatter <- plot1 + plot2
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\featurescatter.jpeg", width=1000, height=350)
feature_scatter
dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)