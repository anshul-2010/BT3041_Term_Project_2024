pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
## Normalizing layer: counts
## Performing log-normalization
## 0%   10   20   30   40   50   60   70   80   90   100%
## [----|----|----|----|----|----|----|----|----|----|
## **************************************************|

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
## Finding variable features for layer counts
## Calculating gene variances
## 0%   10   20   30   40   50   60   70   80   90   100%
## [----|----|----|----|----|----|----|----|----|----|
## **************************************************|
## Calculating feature variances of standardized and clipped values
## 0%   10   20   30   40   50   60   70   80   90   100%
## [----|----|----|----|----|----|----|----|----|----|
## **************************************************|

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
feature_plot <- plot1 + plot2
jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\variablefeature.jpeg",  width=1000, height=350)
feature_plot
dev.off()

all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")