# Load required libraries
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
options(Seurat.object.assay.version = "v5")
# reticulate::use_python("/Users/suoqinjin/anaconda3/bin/python", required=T) 

# R version 4.3.1
# Seurat version 5.0.1

## Data Input
# Subset seurat object by treatment
seurat_ra <- subset(merged_data, subset = Treatment == 'RA')
seurat_sle <- subset(merged_data, subset = Treatment == 'SLE')
seurat_hc <- subset(merged_data, subset = Treatment == 'HC')

# Create CellChat objects from seurat objects
cc_ra <- createCellChat(object = seurat_ra, group.by = "ident", assay = "RNA")
cc_sle <- createCellChat(object = seurat_sle, group.by = "ident", assay = "RNA")
cc_hc <- createCellChat(object = seurat_hc, group.by = "ident", assay = "RNA")

## Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# Only uses the Secreted Signaling from CellChatDB v1
CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))

# set the used database in the object
cc_ra@DB <- CellChatDB.use
cc_sle@DB <- CellChatDB.use
cc_hc@DB <- CellChatDB.use


## Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cc_ra <- subsetData(cc_ra) # This step is necessary even if using the whole database
cc_sle <- subsetData(cc_sle) # This step is necessary even if using the whole database
cc_hc <- subsetData(cc_hc) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cc_ra <- identifyOverExpressedGenes(cc_ra)
cc_ra <- identifyOverExpressedInteractions(cc_ra)
cc_sle <- identifyOverExpressedGenes(cc_sle)
cc_sle <- identifyOverExpressedInteractions(cc_sle)
cc_hc <- identifyOverExpressedGenes(cc_hc)
cc_hc <- identifyOverExpressedInteractions(cc_hc)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cc_ra <- projectData(cc_ra, PPI.human)
cc_sle <- projectData(cc_sle, PPI.human)
cc_hc <- projectData(cc_hc, PPI.human)

### Inference of cell-cell communication network
## Compute the communication probability and infer cellular communication network
cc_ra <- computeCommunProb(cc_ra, type = "triMean", raw.use = "FALSE")
cc_ra <- filterCommunication(cc_ra, min.cells = 10)

cc_sle <- computeCommunProb(cc_sle, type = "triMean", raw.use = "FALSE")
cc_sle <- filterCommunication(cc_sle, min.cells = 10)

cc_hc <- computeCommunProb(cc_hc, type = "triMean", raw.use = "FALSE")
cc_hc <- filterCommunication(cc_hc, min.cells = 10)

## Infer the cell-cell communication at a signaling pathway level
cc_ra <- computeCommunProbPathway(cc_ra)
cc_sle <- computeCommunProbPathway(cc_sle)
cc_hc <- computeCommunProbPathway(cc_hc)

## Calculate the aggregated cell-cell communication network 
cc_ra <- aggregateNet(cc_ra)
cc_sle <- aggregateNet(cc_sle)
cc_hc <- aggregateNet(cc_hc)


### Visualization of cell-cell communication network
### Mif Pathway

# Set pathway
pathways.mif <- c("MIF")

## Visualize signaling strength

## SLE
pairLR.MIF <- extractEnrichedLR(cc_sle, signaling = pathways.mif, geneLR.return = FALSE)
LR.show <- pairLR.MIF[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cc_sle, signaling = pathways.mif,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cc_sle, signaling = pathways.mif, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cc_sle, signaling = pathways.mif, pairLR.use = LR.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cc_sle, signaling = pathways.mif, color.heatmap = "Reds")

## HC
pairLR.MIF <- extractEnrichedLR(cc_hc, signaling = pathways.mif, geneLR.return = FALSE)
LR.show <- pairLR.MIF[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cc_hc, signaling = pathways.mif,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cc_hc, signaling = pathways.mif, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cc_hc, signaling = pathways.mif, pairLR.use = LR.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cc_hc, signaling = pathways.mif, color.heatmap = "Reds")

## Plot the signaling gene expression distribution using violin/dot plot

## SLE
plotGeneExpression(cc_sle, signaling = "MIF", enriched.only = TRUE, type = "violin")

## HC
plotGeneExpression(cc_hc, signaling = "MIF", enriched.only = TRUE, type = "violin")


### Galectin

# Set pathway
pathways.gal <- c("GALECTIN")

## Visualize signaling strength

## RA
pairLR.GALECTIN <- extractEnrichedLR(cc_ra, signaling = pathways.gal, geneLR.return = FALSE)
LR.show <- pairLR.GALECTIN[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cc_ra, signaling = pathways.gal,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cc_ra, signaling = pathways.gal, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cc_ra, signaling = pathways.gal, pairLR.use = LR.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cc_ra, signaling = pathways.gal, color.heatmap = "Reds")

## HC
pairLR.GALECTIN <- extractEnrichedLR(cc_hc, signaling = pathways.gal, geneLR.return = FALSE)
LR.show <- pairLR.GALECTIN[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cc_hc, signaling = pathways.gal,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cc_hc, signaling = pathways.gal, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cc_hc, signaling = pathways.gal, pairLR.use = LR.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cc_hc, signaling = pathways.gal, color.heatmap = "Reds")

## Plot the signaling gene expression distribution using violin/dot plot

## RA
plotGeneExpression(cc_ra, signaling = "GALECTIN", enriched.only = TRUE, type = "violin")

## HC
plotGeneExpression(cc_hc, signaling = "GALECTIN", enriched.only = TRUE, type = "violin")
