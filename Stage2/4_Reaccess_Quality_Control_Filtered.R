# RE_ASSESS QC METRICS
View(filtered_seurat@meta.data)

# Add number of genes per UMI for each cell to metadata
filtered_seurat$log10GenesPerUMI <- log10(filtered_seurat$nFeature_RNA) / log10(filtered_seurat$nCount_RNA)

# Compute percent mito ratio
filtered_seurat$mitoRatio <- PercentageFeatureSet(object = filtered_seurat, pattern = "^MT-")
filtered_seurat$mitoRatio <- filtered_seurat@meta.data$mitoRatio / 100

# Create metadata frame
metadata <- filtered_seurat@meta.data
metadata$cells <- rownames(metadata)

metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"

metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,nUMI = nCount_RNA,nGene = nFeature_RNA)

metadata$mitoRatio <- filtered_seurat@meta.data$mitoRatio

View(metadata)

filtered_seurat@meta.data <- metadata
save(filtered_seurat, file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\GSE262446_RAW\\filtered_filtered_seurat.RData")

jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\cell_counts_filtered.jpeg", width=500)
metadata %>% 
  	ggplot(aes(x=sample, fill=sample)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")
dev.off()

jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\umi_counts_per_cell_filtered.jpeg", width=500)
metadata %>% 
  	ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500)
dev.off()

jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\genes_per_cell_filtered.jpeg", width=500)
metadata %>% 
  	ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)
dev.off()

jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\complexity_filtered.jpeg", width=500)
metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)
dev.off()

jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\mito_counts_ratio_filtered.jpeg", width=500)
metadata %>% 
  	ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)
dev.off()

jpeg(file="C:\\Users\\Dell\\Desktop\\Courses\\Sem_VI\\BT3041\\Project\\All_Plots\\joint_filtering_effect_filtered.jpeg", width=500)
metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~sample)
dev.off()