library(dplyr)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(ggplot2)
library(patchwork)
library(hdf5r)
library(tidyr)


# R version 4.3.1
# Seurat version 5.0.1

# Import data and create seurat objects
RA.sample <- Read10X_h5("GSM4819747_RA_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE) # check
RA <-CreateSeuratObject(RA.sample,  project = "RA", min.cells = 3)

SLE.1 <- Read10X("GSM4954811_SLE_1")
SLE1 <- CreateSeuratObject(SLE.1, project = "SLE_1", min.cells = 3)

SLE.2 <- Read10X("GSM4954812_SLE-2")
SLE2 <- CreateSeuratObject(SLE.2, project = "SLE_2", min.cells = 3)

SLE.3 <- Read10X("GSM4217720_scRNA_SLE3")
SLE3 <- CreateSeuratObject(SLE.3, project = "SLE_3", min.cells = 3)

HC.1 <- Read10X("GSM4954813_C-1")
HC1 <- CreateSeuratObject(HC.1, project = "HC_1", min.cells = 3)

HC.2 <- Read10X("GSM4954813_C-2")
HC2 <- CreateSeuratObject(HC.2, project = "HC_2", min.cells = 3)

HC.3 <- Read10X("GSM4954813_C-3")
HC3 <- CreateSeuratObject(HC.3, project = "HC_3", min.cells = 3)


# Merge all datasets
merged_data <- merge(RA, y = c(SLE1, SLE2, SLE3, HC1), 
                     add.cell.ids = c("RA_1","SLE_1", "SLE_2", "SLE_3", "HC_1"), 
                     project = "autoimmune")
merged_data
table(merged_data$orig.ident)

merged_data$sample <- rownames(merged_data@meta.data)

# split sample column
merged_data@meta.data <- separate(merged_data@meta.data, col = 'sample', 
                                  into = c("Treatment","Patient", "Barcode"), sep = "_")

# Add column to metadata table looking at mitochondrial percentage
mito_genes <- c("MT-ND1", "MT-ND2","MT-CO1","MT-CO2","MT-ATP8","MT-ATP6","MT-CO3", 
                "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB")
merged_data[["percent.mt"]] <- PercentageFeatureSet(object=merged_data, features=mito_genes, assay="RNA")

# Visualize QC metrics as a violin plot
VlnPlot(merged_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(merged_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter the dataset
merged_data <- subset(merged_data, subset = nFeature_RNA > 200 & percent.mt < 10)
merged_data

# run standard anlaysis workflow
merged_data <- NormalizeData(merged_data)
merged_data <- FindVariableFeatures(merged_data)
merged_data <- ScaleData(merged_data)
merged_data <- RunPCA(merged_data)

# Check to see how many PCs should be used
ElbowPlot(merged_data, ndims = 50)

# Perform integration
merged_data <- IntegrateLayers(object = merged_data, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE, features = VariableFeatures(merged_data))

# re-join layers after integration
merged_data[["RNA"]] <- JoinLayers(merged_data[["RNA"]])

# Create clusters
merged_data <- FindNeighbors(merged_data, reduction = "integrated.cca", dims = 1:35)
merged_data <- FindClusters(merged_data, resolution = 0.6)
merged_data <- RunTSNE(merged_data, dims = 1:35, reduction = "integrated.cca")

# Visualization of clusters and marker genes
DimPlot(merged_data, reduction = "tsne", group.by = "Treatment")
DimPlot(merged_data, reduction = "tsne", label = T)
FeaturePlot(merged_data, features = c("MS4A1", "CD8A", "IL7R", "CD14"), max.cutoff = 6, cols = c("grey",                                                                                              "red"), reduction = "tsne")
FeaturePlot(merged_data, features = c("FCGR3A", "NKG7", "FCER1A", "PPBP"), max.cutoff = 6, cols = c("grey",
                                                                                                 "red"), reduction = "tsne")

# Annotate clusters
merged_data <- RenameIdents(merged_data, `0` = "CD14+ mono", `1` = "NK", `2` = "CD4+ T", `4` = "CD8+ T",
                            `3` = "mixed", `5` = "CD4+ T", `6` = "NK", `7` = "B",'8' ="B", `9` = "FCGR3A+ mono", `10` = "CD4+ T",
                            `11` = "CD4+ T", `12` = "mixed", `13` = "NK", `14` = "CD4+ T", `15` = "platelet", 
                            `16` = "mixed",`17` = "DC", `18`='mixed', `19`="mixed", `20`="mixed")
DimPlot(merged_data, label = TRUE)



# create a meta.data column for celltypes and celltypes_treatment
merged_data$celltype.treatment <- paste(Idents(merged_data), merged_data$Treatment, sep = "_")
merged_data$celltypes <- Idents(merged_data)
Idents(merged_data) <- "celltypes"

# create a table that contains counts for each cell type. These are used to create
# a stacked bar chart that contains proportions for each cluster
counts_by_patient <- table(merged_data@meta.data$celltypes, merged_data@meta.data$orig.ident)
write.csv(counts_by_patient, "proportions.csv", row.names = TRUE)

# Identify the top 3 genes that define each cluster
nk.markers <- FindConservedMarkers(merged_data, ident.1 = "NK", grouping.var = "Treatment", verbose = FALSE)
head(nk.markers, n =3)
fcgr3a.markers <- FindConservedMarkers(merged_data, ident.1 = "FCGR3A+ mono", grouping.var = "Treatment", verbose = FALSE)
head(fcgr3a.markers, n =3)
b.markers <- FindConservedMarkers(merged_data, ident.1 = "B", grouping.var = "Treatment", verbose = FALSE)
head(b.markers, n =3)
cd8.markers <- FindConservedMarkers(merged_data, ident.1 = "CD8+ T", grouping.var = "Treatment", verbose = FALSE)
head(cd8.markers, n =3)
cd4.markers <- FindConservedMarkers(merged_data, ident.1 = "CD4+ T", grouping.var = "Treatment", verbose = FALSE)
head(cd4.markers, n =3)
platelet.markers <- FindConservedMarkers(merged_data, ident.1 = "platelet", grouping.var = "Treatment", verbose = FALSE)
head(platelet.markers, n =3)
dc.markers <- FindConservedMarkers(merged_data, ident.1 = "DC", grouping.var = "Treatment", verbose = FALSE)
head(dc.markers, n =3)
cd14.markers <- FindConservedMarkers(merged_data, ident.1 = "CD14+ mono", grouping.var = "Treatment", verbose = FALSE)
head(cd14.markers, n =3)

# Create heatmap that identifies the top 10 genes that define each cluster
merged.markers <- FindAllMarkers(merged_data, only.pos = TRUE)
merged.markers %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 1)
merged.markers %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 1) %>%
        slice_head(n = 10) %>%
        ungroup() -> top10
groups_of_interest <- subset(merged_data, idents = c("CD14+ mono", "CD4+ T", "CD8+ T", "platelet", "DC", "NK", "FCGR3A+ mono", "B"))
groups_of_interest$celltypes_minus_mixed <- Idents(groups_of_interest)
map <-DoHeatmap(groups_of_interest, features = top10$gene, group.by = "celltypes_minus_mixed", angle = 45, hjust = 0, size = 3.0) + 
        scale_fill_gradient(low = "lightblue", high = "red", na.value = NA) + theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5))


# Identify differentialy expressed genes between cell types (SLE vs HC) 
# ribosomal and ATP genes were filtered out
B.SLE.HC <- FindMarkers(merged_data, ident.1 = "B_SLE", ident.2 = "B_HC", verbose = FALSE, 
                        min.pct = 0.25, logfc.threshold = 0.25,fdr.threshold = 0.05,
                        only.pos = T)
ribosomal_genes <- grep("^RPL|^RPS", rownames(B.SLE.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(B.SLE.HC), value = TRUE)
B.SLE.HC.filtered_results <- B.SLE.HC[!(rownames(B.SLE.HC) %in% ribosomal_genes), ]
B.SLE.HC.filtered_results <- B.SLE.HC.filtered_results[!(rownames(B.SLE.HC.filtered_results) %in% ATP_genes), ]
write.csv(B.SLE.HC.filtered_results, "B.SLE.HC.test.csv", row.names=TRUE)

CD8.SLE.HC <- FindMarkers(merged_data, ident.1 = "CD8+ T_SLE", ident.2 = "CD8+ T_HC", 
                          verbose = FALSE, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE, fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(CD8.SLE.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(CD8.SLE.HC), value = TRUE)
CD8.SLE.HC.filtered_results <- CD8.SLE.HC[!(rownames(CD8.SLE.HC) %in% ribosomal_genes), ]
CD8.SLE.HC.filtered_results <- CD8.SLE.HC.filtered_results[!(rownames(CD8.SLE.HC.filtered_results) %in% ATP_genes), ]
write.csv(CD8.SLE.HC.filtered_results, "CD8.SLE.HC.csv", row.names=TRUE)

CD14.SLE.HC <- FindMarkers(merged_data, ident.1 = "CD14+ mono_SLE", ident.2 = "CD14+ mono_HC", 
                           verbose = FALSE, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE,fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(CD14.SLE.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(CD14.SLE.HC), value = TRUE)
CD14.SLE.HC.filtered_results <- CD14.SLE.HC[!(rownames(CD14.SLE.HC) %in% ribosomal_genes), ]
CD14.SLE.HC.filtered_results <- CD14.SLE.HC.filtered_results[!(rownames(CD14.SLE.HC.filtered_results) %in% ATP_genes), ]
write.csv(CD14.SLE.HC.filtered_results, "CD14.SLE.HC.csv", row.names=TRUE)

CD4.SLE.HC <- FindMarkers(merged_data, ident.1 = "CD4+ T_SLE", ident.2 = "CD4+ T_HC", 
                          verbose = FALSE, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE,fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(CD4.SLE.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(CD4.SLE.HC), value = TRUE)
CD4.SLE.HC.filtered_results <- CD4.SLE.HC [!(rownames(CD4.SLE.HC ) %in% ribosomal_genes), ]
CD4.SLE.HC.filtered_results <-CD4.SLE.HC.filtered_results[!(rownames(CD4.SLE.HC.filtered_results) %in% ATP_genes), ]
write.csv(CD4.SLE.HC.filtered_results, "CD4.SLE.HC.csv", row.names=TRUE)

DC.SLE.HC <- FindMarkers(merged_data, ident.1 = "DC_SLE", ident.2 = "DC_HC", 
                         verbose = FALSE, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE,fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(DC.SLE.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(DC.SLE.HC), value = TRUE)
DC.SLE.HC.filtered_results <- DC.SLE.HC[!(rownames(DC.SLE.HC) %in% ribosomal_genes), ]
DC.SLE.HC.filtered_results <- DC.SLE.HC.filtered_results[!(rownames(DC.SLE.HC.filtered_results) %in% ATP_genes), ]
write.csv(DC.SLE.HC.filtered_results, "DC.SLE.HC.csv", row.names=TRUE)

FCG.SLE.HC <- FindMarkers(merged_data, ident.1 = "FCGR3A+ mono_SLE", ident.2 = "FCGR3A+ mono_HC", 
                          verbose = FALSE, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE,fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(FCG.SLE.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(FCG.SLE.HC), value = TRUE)
FCG.SLE.HC.filtered_results <- FCG.SLE.HC[!(rownames(FCG.SLE.HC) %in% ribosomal_genes), ]
FCG.SLE.HC.filtered_results <- FCG.SLE.HC.filtered_results[!(rownames(FCG.SLE.HC.filtered_results) %in% ATP_genes), ]
write.csv(FCG.SLE.HC.filtered_results , "FCG.SLE.HC.csv", row.names=TRUE)

NK.SLE.HC <- FindMarkers(merged_data, ident.1 = "NK_SLE", ident.2 = "NK_HC", 
                         verbose = FALSE, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE,fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(NK.SLE.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(NK.SLE.HC), value = TRUE)
NK.SLE.HC.filtered_results <- NK.SLE.HC[!(rownames(NK.SLE.HC) %in% ribosomal_genes), ]
NK.SLE.HC.filtered_results <- NK.SLE.HC.filtered_results[!(rownames(NK.SLE.HC.filtered_results) %in% ATP_genes), ]
write.csv(NK.SLE.HC.filtered_results, "NK.SLE.HC.csv", row.names=TRUE)

platelet.SLE.HC <- FindMarkers(merged_data, ident.1 = "platelet_SLE", ident.2 = "platelet_HC", 
                               verbose = FALSE, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE,fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(platelet.SLE.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(platelet.SLE.HC), value = TRUE)
platelet.SLE.HC.filtered_results <- platelet.SLE.HC[!(rownames(platelet.SLE.HC) %in% ribosomal_genes), ]
platelet.SLE.HC.filtered_results <- platelet.SLE.HC.filtered_results[!(rownames(platelet.SLE.HC.filtered_results) %in% ATP_genes), ]
write.csv(platelet.SLE.HC.filtered_results, "platelet.SLE.HC.csv", row.names=TRUE)


# Identify differentialy expressed genes between each cell type (RA vs HC) 
# ribosomal and ATP genes were filtered out
B.RA.HC <- FindMarkers(merged_data, ident.1 = "B_RA", ident.2 = "B_HC", verbose = FALSE, 
                        min.pct = 0.25, logfc.threshold = 1.0,fdr.threshold = 0.05,
                        only.pos = T)
ribosomal_genes <- grep("^RPL|^RPS", rownames(B.RA.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(B.RA.HC), value = TRUE)
B.RA.HC.filtered_results <- B.RA.HC[!(rownames(B.RA.HC) %in% ribosomal_genes), ]
B.RA.HC.filtered_results <- B.RA.HC.filtered_results[!(rownames(B.RA.HC.filtered_results) %in% ATP_genes), ]
write.csv(B.RA.HC.filtered_results, "B.RA.HC.csv", row.names=TRUE)

CD8.RA.HC <- FindMarkers(merged_data, ident.1 = "CD8+ T_RA", ident.2 = "CD8+ T_HC", 
                          verbose = FALSE, min.pct = 0.25, logfc.threshold = 1.0, only.pos = TRUE, fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(CD8.RA.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(CD8.RA.HC), value = TRUE)
CD8.RA.HC.filtered_results <- CD8.RA.HC[!(rownames(CD8.RA.HC) %in% ribosomal_genes), ]
CD8.RA.HC.filtered_results <- CD8.RA.HC.filtered_results[!(rownames(CD8.RA.HC.filtered_results) %in% ATP_genes), ]
write.csv(CD8.RA.HC.filtered_results, "CD8.RA.HC.csv", row.names=TRUE)

CD14.RA.HC <- FindMarkers(merged_data, ident.1 = "CD14+ mono_RA", ident.2 = "CD14+ mono_HC", 
                           verbose = FALSE, min.pct = 0.25, logfc.threshold = 1.0, only.pos = TRUE,fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(CD14.RA.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(CD14.RA.HC), value = TRUE)
CD14.RA.HC.filtered_results <- CD14.RA.HC[!(rownames(CD14.RA.HC) %in% ribosomal_genes), ]
CD14.RA.HC.filtered_results <- CD14.RA.HC.filtered_results[!(rownames(CD14.RA.HC.filtered_results) %in% ATP_genes), ]
write.csv(CD14.RA.HC.filtered_results, "CD14.RA.HC.csv", row.names=TRUE)

CD4.RA.HC <- FindMarkers(merged_data, ident.1 = "CD4+ T_RA", ident.2 = "CD4+ T_HC", 
                          verbose = FALSE, min.pct = 0.25, logfc.threshold = 1.0, only.pos = TRUE,fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(CD4.RA.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(CD4.RA.HC), value = TRUE)
CD4.RA.HC.filtered_results <- CD4.RA.HC[!(rownames(CD4.RA.HC) %in% ribosomal_genes), ]
CD4.RA.HC.filtered_results <- CD4.RA.HC.filtered_results[!(rownames(CD4.RA.HC.filtered_results) %in% ATP_genes), ]
write.csv(CD4.RA.HC.filtered_results, "CD4.RA.HC.csv", row.names=TRUE)

DC.RA.HC <- FindMarkers(merged_data, ident.1 = "DC_RA", ident.2 = "DC_HC", 
                         verbose = FALSE, min.pct = 0.25, logfc.threshold = 1.0, only.pos = TRUE,fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(DC.RA.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(DC.RA.HC), value = TRUE)
DC.RA.HC.filtered_results <- DC.RA.HC[!(rownames(DC.RA.HC) %in% ribosomal_genes), ]
DC.RA.HC.filtered_results <- DC.RA.HC.filtered_results[!(rownames(DC.RA.HC.filtered_results) %in% ATP_genes), ]
write.csv(DC.RA.HC.filtered_results, "DC.RA.HC.csv", row.names=TRUE)

FCG.RA.HC <- FindMarkers(merged_data, ident.1 = "FCGR3A+ mono_RA", ident.2 = "FCGR3A+ mono_HC", 
                          verbose = FALSE, min.pct = 0.25, logfc.threshold = 1.0, only.pos = TRUE,fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(FCG.RA.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(FCG.RA.HC), value = TRUE)
FCG.RA.HC.filtered_results <- FCG.RA.HC[!(rownames(FCG.RA.HC) %in% ribosomal_genes), ]
FCG.RA.HC.filtered_results <- FCG.RA.HC.filtered_results[!(rownames(FCG.RA.HC.filtered_results) %in% ATP_genes), ]
write.csv(FCG.RA.HC.filtered_results, "FCG.RA.HC.csv", row.names=TRUE)

NK.RA.HC <- FindMarkers(merged_data, ident.1 = "NK_RA", ident.2 = "NK_HC", 
                         verbose = FALSE, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE,fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(NK.RA.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(NK.RA.HC), value = TRUE)
NK.RA.HC.filtered_results <- NK.RA.HC[!(rownames(NK.RA.HC) %in% ribosomal_genes), ]
NK.RA.HC.filtered_results <- NK.RA.HC.filtered_results[!(rownames(NK.RA.HC.filtered_results) %in% ATP_genes), ]
write.csv(NK.RA.HC.filtered_results, "NK.RA.HC.csv", row.names=TRUE)

platelet.RA.HC <- FindMarkers(merged_data, ident.1 = "platelet_RA", ident.2 = "platelet_HC", 
                               verbose = FALSE, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE,fdr.threshold = 0.05)
ribosomal_genes <- grep("^RPL|^RPS", rownames(platelet.RA.HC), value = TRUE)
ATP_genes <- grep("^ATP", rownames(platelet.RA.HC), value = TRUE)
platelet.RA.HC.filtered_results <- platelet.RA.HC[!(rownames(platelet.RA.HC) %in% ribosomal_genes), ]
platelet.RA.HC.filtered_results <- platelet.RA.HC.filtered_results[!(rownames(platelet.RA.HC.filtered_results) %in% ATP_genes), ]
write.csv(platelet.RA.HC.filtered_results, "platelet.RA.HC.csv", row.names=TRUE)
