rm(list=ls())

setwd("/Users/svalpione/Dropbox (The University of Manchester)/projects/Single cell TME/GSE100501_RAW")
library(Seurat)
library(calibrate)
library(data.table)


#saveRDS(immune.combined, file = "immune.combined.rds")


install.packages('devtools')
devtools::install_version(package = 'Seurat', version = package_version('3.0.1'))

install.packages('devtools')
devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))


library(Matrix)
library(fieldRS)
y
#comparison between tot cells and Stim cells

reap_rna_2 <- read.table("GSM2685238_mRNA_2_PBMCs_matrix.txt.gz", row.names = 1)
reap_rna_3 <- read.table("GSM2685239_mRNA_3_PBMCs_matrix.txt.gz", row.names = 1)

reap_protein_2 <- read.table("GSM2685243_protein_2_PBMCs_matrix.txt.gz", row.names = 1)
reap_protein_3 <- read.table("GSM2685244_protein_3_PBMCs_matrix.txt.gz", row.names = 1)


setwd("/Users/svalpione/Dropbox (The University of Manchester)/projects/Single cell TME/GSE100501_RAW/results26May19")

reap_rna_combined <- Matrix(as.matrix(cbind(reap_rna_2, reap_rna_3)), sparse = TRUE)
reap_protein_combined <- Matrix(as.matrix(cbind(reap_protein_2, reap_protein_3)), sparse = TRUE)

rownames(reap_protein_combined) <- make.unique(paste0("REAP_", sapply(rownames(reap_protein_combined), ExtractField)))

reap_protein_combined <- reap_protein_combined[setdiff(rownames(reap_protein_combined), c("REAP_Control_Rat_IgG1","REAP_Blank", "REAP_Control", "REAP_Control.1", "REAP_Control_Mouse_IgG2b_CATAAAGG","REAP_Mouse_IgG1_TGTGTATA","REAP_Control_Rat_IgG1_TCACGGTA")), ]

library(Seurat)
head(reap_protein_combined)

reap <- CreateSeuratObject(raw.data = reap_rna_combined, project = "IMMUNE_CTRL")

reap@meta.data$Stim <- "CTRL"

reap <- NormalizeData(reap)
reap <- FindVariableFeatures(reap)
reap <- ScaleData(reap, display.progress = F)

reap<- RunPCA(reap, verbose=FALSE)
reap <- RunTSNE(reap, dims = 1:25, method = "FIt-SNE")


reap[["REAP"]]<-CreateAssayObject(counts=reap_protein_combined)
reap <- NormalizeData(reap, assay.type = "REAP", normalization.method = "CLR", display.progress = FALSE)
reap <- ScaleData(reap, assay.type = "REAP", display.progress = TRUE)


FeaturePlot(reap, features = c("CD3"), min.cutoff = "q05", max.cutoff = "q95", ncol = 4, pt.size = 0.5)


reap <- SetAssayData(reap, assay.type = "REAP", slot = "counts", new.data = reap_protein_combined[, reap@counts])



reap_REAP<- RunPCA(reap, pc.genes = rownames(reap_protein_combined), assay.type = "REAP", pcs.print = 0)
PCAPlot(reap_REAP, pt.size = 0.5)
adt.data <- GetAssayData(reap_REAP, assay.type = "REAP", slot = "data")
reap_REAP <- StashIdent(reap_REAP, "rnaClusterID")
reap_REAP <- RunTSNE(reap_REAP, dims.use = 1:13)
TSNEPlot(object = reap_REAP)

FeaturePlot(reap_REAP, features.plot = c("REAP_CD3_AGGATCGA", "REAP_CD8a_ACCCGCAC", "REAP_CD4_CACGATTC", "REAP_CD45RA_GTGATAGT", "REAP_CD45RO_TGATATCG", "REAP_CD27_GCTGTGTA",REAP_CD197_ACGCTTGG), min.cutoff = "q05", max.cutoff = "q95", nCol = 4, cols.use = c("lightgrey", "blue"), pt.size = 0.5)



summary(FetchData(reap_REAP, vars.all='REAP_CD45RO'))
summary(FetchData(reap_REAP, vars.all='REAP_CD45RA'))
summary(FetchData(reap_REAP, vars.all='REAP_CD27'))
summary(FetchData(reap_REAP, vars.all='REAP_CD197'))
summary(FetchData(reap_REAP, vars.all='REAP_CD8a'))

REAP_CD45RO.cutoff <- 0.66
REAP_CD45RA.cutoff <- 0.63
REAP_CD27.cutoff <- 0
REAP_CD197.cutoff <- 0
REAP_CD8a.cutoff <- 1

length(which(FetchData(reap_REAP, vars.all ='REAP_CD45RO') >= REAP_CD45RO.cutoff))
length(which(FetchData(reap_REAP, vars.all='REAP_CD27') <= REAP_CD27.cutoff))
length(which(FetchData(reap_REAP, vars.all='REAP_CD197') <= REAP_CD197.cutoff))
length(which(FetchData(reap_REAP, vars.all='REAP_CD8a') >= REAP_CD8a.cutoff))
length(which(FetchData(reap_REAP, vars.all='REAP_CD45RA') <= REAP_CD45RA.cutoff))

length(which(FetchData(reap_REAP, vars.all='REAP_CD45RO') >= REAP_CD45RO.cutoff & FetchData(reap_REAP, vars.all='REAP_CD45RA') <= REAP_CD45RA.cutoff & FetchData(reap_REAP, vars.all='REAP_CD27') <= REAP_CD27.cutoff & FetchData(reap_REAP, vars.all='REAP_CD197') <= REAP_CD197.cutoff & FetchData(reap_REAP, vars.all='REAP_CD8a') >= REAP_CD8a.cutoff ))

C<-(which(FetchData(reap_REAP, vars.all='REAP_CD45RO') >= REAP_CD45RO.cutoff & FetchData(reap_REAP, vars.all='REAP_CD45RA') <= REAP_CD45RA.cutoff & FetchData(reap_REAP, vars.all='REAP_CD27') <= REAP_CD27.cutoff & FetchData(reap_REAP, vars.all='REAP_CD197') <= REAP_CD197.cutoff & FetchData(reap_REAP, vars.all='REAP_CD8a') >= REAP_CD8a.cutoff))

tmp.t <- SubsetData(object = reap_REAP, cells.use = reap_REAP@cell.names[C])
tmp.t
TMP.t<-c(tmp.t@cell.names)

#reap_REAP <- SetIdent(object = reap_REAP, cells.use = TMP.t, ident.use = "Tcruk.t")


quartz.save("REAP_ridge_plot_stot_Tcruk.pdf", type="pdf")


reap_REAP <- FindVariableGenes(reap_REAP, do.plot = F)
stim_REAP <- FindVariableGenes(stim_REAP, do.plot = F)
g.1 <- head(rownames(reap_REAP@hvg.info), 1000)
g.2 <- head(rownames(stim_REAP@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(reap_REAP@scale.data))
genes.use <- intersect(genes.use, rownames(stim_REAP@scale.data))

immune.combined <- RunCCA(reap_REAP, stim_REAP, num.cc = 30)


# visualize results of CCA plot CC1 versus CC2 and look at a violin plot

p1 <- DimPlot(object = immune.combined, reduction.use = "cca", group.by = "Stim", pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "CC1", group.by = "Stim", do.return = TRUE)
plot_grid(p1, p2)
quartz.save("immune.combined.pdf", type="PDF")



PrintDim(object = immune.combined, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

DimHeatmap(object = immune.combined, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)

immune.combined <- CalcVarExpRatio(object = immune.combined, reduction.type = "pca", grouping.var = "Stim", dims.use = 1:13)

immune.combined.all.save <- immune.combined
immune.combined <- SubsetData(object = immune.combined, subset.name = "var.ratio.pca", accept.low = 0.5)

immune.combined.discard <- SubsetData(object = immune.combined.all.save, subset.name = "var.ratio.pca", accept.high = 0.5)

immune.combined <- AlignSubspace(object = immune.combined, reduction.type = "cca", grouping.var ="Stim", dims.align = 1:13)

p1 <- VlnPlot(object = immune.combined, features.plot = "ACC1", group.by = "Stim", 
    do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "ACC2", group.by = "Stim", 
    do.return = TRUE)
plot_grid(p1, p2)
quartz.save("immune.combined.AAC.pdf", type="PDF")

dev.off()


immune.combined <- RunTSNE(object = immune.combined, reduction.use = "cca.aligned", dims.use = 1:13, do.fast = TRUE)
immune.combined <- FindClusters(object = immune.combined, reduction.type = "cca.aligned", dims.use = 1:13, save.SNN = TRUE)


immune.combined <- SetIdent(object = immune.combined, cells.use = TMP.s, ident.use = "Tcruk.s")
immune.combined <- SetIdent(object = immune.combined, cells.use = TMP.t, ident.use = "Tcruk.t")

immune.combined@dr



p1 <- TSNEPlot(object = immune.combined, group.by = "Stim", do.return = TRUE, pt.size = 0.5)
p2 <- TSNEPlot(object = immune.combined, do.return = TRUE, pt.size = 0.5)
plot_grid(p1, p2)
quartz.save("immune.combined.tot.tSNE_Tcruk.pdf", type="PDF")

immune.combined.markers <- FindMarkers(object = immune.combined, ident.1 = "Tcruk.t",  ident.2 = "Tcruk.s", assay.type = "RNA", logfc.threshold=0.25, only.pos = FALSE)
head(immune.combined.markers)

write.csv(immune.combined.markers, "immune.combined.markers_Tcruk_post_stim.tot.csv")
vol<-read.csv("immune.combined.markers_Tcruk_post_stim.tot.csv", sep=",")
head(vol)

vol$p_val[vol$p_val==0]<-10^-6

GenePlot(immune.combined, gene1="GZMA", gene2="KLRG1")

with(vol, plot(avg_logFC, -log10(p_val), pch=20, main="Tcruk differential expression after stimulation",xlim=c(-2.5,2.4)))

with(subset(vol, p_val_adj <.05 ), points(avg_logFC, -log10(p_val), pch=20, col="red"))
with(subset(vol, abs(avg_logFC)>0.5), points(avg_logFC, -log10(p_val), pch=20, col="orange"))
with(subset(vol, p_val_adj <.05 & abs(avg_logFC)>1), points(avg_logFC, -log10(p_val), pch=20, col="green"))


with(subset(vol, p_val_adj <.05 & abs(avg_logFC)>1), textxy(avg_logFC, -log10(p_val), labs=X, cex=.4))

quartz.save("Tcruk differential expression after stimulation_tot.pdf", type="PDF")


saveRDS(immune.combined, file = "immune.combined.rds")















#differential expression Tcruk vs all the rest for CD69 and EOMES
reap <- CreateSeuratObject(raw.data = reap_rna_combined)
reap <- NormalizeData(reap)
reap <- FindVariableGenes(reap, do.plot = FALSE, y.cutoff = 0.5)
reap <- ScaleData(reap, display.progress = FALSE)
reap <- RunPCA(reap, pcs.print = 0)
reap <- FindClusters(reap, dims.use = 1:13, print.output = FALSE)
reap <- RunTSNE(reap, dims.use = 1:13)
reap.rna.markers <- FindAllMarkers(reap, max.cells.per.ident = 100, logfc.threshold = log(2), only.pos = TRUE, min.diff.pct = 0.3, do.print = F)
TSNEPlot(reap, do.label = TRUE, pt.size = 0.5)
head(reap.rna.markers)

reap <- SetAssayData(reap, assay.type = "REAP", slot = "raw.data", new.data = reap_protein_combined)
reap <- NormalizeData(reap, assay.type = "REAP", normalization.method = "genesCLR", display.progress = FALSE)
reap <- ScaleData(reap, assay.type = "REAP",do.scale=TRUE, display.progress = TRUE)
FeaturePlot(reap, features.plot = c("REAP_CD3", "REAP_CD8", "REAP_CD4", "REAP_CD45RA","REAP_CD45RO", "REAP_CD27","REAP_CD197"), min.cutoff = "q05", max.cutoff = "q95", nCol = 4, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
RidgePlot(reap_REAP, features.plot = c("REAP_CD3", "REAP_CD8", "REAP_CD4", "REAP_CD45RA","REAP_CD45RO", "REAP_CD27","REAP_CD197"), nCol = 2)#Tcruk is cluster 0?
prot.markers <- FindMarkers(reap, ident.1="1", assay.type = "REAP", logfc.threshold = log(1.5))#es for cluster 1
head(prot.markers)
adt.markers <- FindAllMarkers(reap, assay.type = "REAP", only.pos = TRUE, print.bar = F)
DoHeatmap(reap, genes.use = unique(adt.markers$gene), assay.type = "REAP", slim.col.label = TRUE, remove.key = TRUE, group.label.rot = TRUE)

reap_REAP<- RunPCA(reap, pc.genes = rownames(reap_protein_combined), assay.type = "REAP", pcs.print = 0)
PCAPlot(reap_REAP, pt.size = 0.5)
adt.data <- GetAssayData(reap_REAP, assay.type = "REAP", slot = "data")
reap_REAP <- StashIdent(reap_REAP, "rnaClusterID")
reap_REAP <- RunTSNE(reap_REAP, dims.use = 1:13)
reap_REAP <- FindClusters(reap_REAP, dims.use = 1:13, print.output = FALSE, force.recalc=TRUE)
clustering.table <- table(reap_REAP@ident, reap_REAP@meta.data$rnaClusterID)
tsne_rnaClusters <- TSNEPlot(reap_REAP, do.return = TRUE, group.by = "rnaClusterID", pt.size = 0.5)
tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Clustering based on scRNA-seq") + theme(plot.title = element_text(hjust = 0.5))
tsne_adtClusters <- TSNEPlot(reap_REAP, do.return = TRUE, pt.size = 0.5)
tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on ADT signal") + theme(plot.title = element_text(hjust = 0.5))
tsne_adtClusters <- TSNEPlot(reap_REAP, do.return = TRUE, pt.size = 0.5)
plot_grid(tsne_rnaClusters, tsne_adtClusters, ncol = 2)
FeaturePlot(reap_REAP, features.plot = c("REAP_CD3", "REAP_CD8", "REAP_CD4", "REAP_CD45RA","REAP_CD45RO", "REAP_CD27","REAP_CD197"), min.cutoff = "q05", max.cutoff = "q95", nCol = 4, cols.use = c("lightgrey", "blue"), pt.size = 0.5)
quartz.save("REAP_feature_plot_tot_Tcruk.pdf", type="pdf")
head(reap@meta.data)


REAP_CD45RO.cutoff <- 0.66
REAP_CD45RA.cutoff <- 0.63
REAP_CD27.cutoff <- 0
REAP_CD197.cutoff <- 0
REAP_CD8a.cutoff <- 1

length(which(FetchData(reap_REAP, vars.all ='REAP_CD45RO') >= REAP_CD45RO.cutoff))
length(which(FetchData(reap_REAP, vars.all='REAP_CD27') <= REAP_CD27.cutoff))
length(which(FetchData(reap_REAP, vars.all='REAP_CD197') <= REAP_CD197.cutoff))
length(which(FetchData(reap_REAP, vars.all='REAP_CD8a') >= REAP_CD8a.cutoff))
length(which(FetchData(reap_REAP, vars.all='REAP_CD45RA') <= REAP_CD45RA.cutoff))

length(which(FetchData(reap_REAP, vars.all='REAP_CD45RO') >= REAP_CD45RO.cutoff & FetchData(reap_REAP, vars.all='REAP_CD45RA') <= REAP_CD45RA.cutoff & FetchData(reap_REAP, vars.all='REAP_CD27') <= REAP_CD27.cutoff & FetchData(reap_REAP, vars.all='REAP_CD197') <= REAP_CD197.cutoff & FetchData(reap_REAP, vars.all='REAP_CD8a') >= REAP_CD8a.cutoff ))

C<-(which(FetchData(reap_REAP, vars.all='REAP_CD45RO') >= REAP_CD45RO.cutoff & FetchData(reap_REAP, vars.all='REAP_CD45RA') <= REAP_CD45RA.cutoff & FetchData(reap_REAP, vars.all='REAP_CD27') <= REAP_CD27.cutoff & FetchData(reap_REAP, vars.all='REAP_CD197') <= REAP_CD197.cutoff & FetchData(reap_REAP, vars.all='REAP_CD8a') >= REAP_CD8a.cutoff))

tmp.t <- SubsetData(object = reap_REAP, cells.use = reap_REAP@cell.names[C])
tmp.t
TMP.t<-c(tmp.t@cell.names)

reap_REAP <- SetIdent(object = reap_REAP, cells.use = TMP.t, ident.use = "Tcruk.t")

#stop, then select the cells...
tsne_adtClusters <- TSNEPlot(reap_REAP, do.return = TRUE, pt.size = 0.5)
tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on ADT signal") + theme(plot.title = element_text(hjust = 0.5))
tsne_adtClusters <- TSNEPlot(reap_REAP, do.return = TRUE, pt.size = 0.5)
plot_grid(tsne_rnaClusters, tsne_adtClusters, ncol = 2)
quartz.save("clusters_RNAvsADT_with Tcruk_tot.pdf", type="PDF")

newcells.markers <- FindMarkers(object = reap_REAP, ident.1 = "Tcruk.t",  assay.type = "RNA", logfc.threshold=0.25, only.pos = TRUE)

head(newcells.markers)
newcells.markers <- FindAllMarkers(object = reap_REAP)

head(x = newcells.markers)


write.csv(newcells.markers, "reap_newcells.markers_Tcruk.csv")

RidgePlot(reap_REAP, features.plot = c("REAP_CD3", "REAP_CD8a", "REAP_CD4", "REAP_CD45RA","REAP_CD45RO", "REAP_CD27","REAP_CD197"), nCol = 2)
quartz.save("reap_tot_ridgePlot_Tcruk.pdf", type="PDF")




reap_REAP <- FindVariableGenes(reap_REAP, do.plot = F)



TSNEPlot(object = reap_REAP, do.return = TRUE, pt.size = 0.5)

immune.markers <- FindMarkers(object = reap_REAP, ident.1 = "Tcruk.t",  ident.2 = "0", assay.type = "RNA", logfc.threshold=0.25, only.pos = FALSE)
head(immune.markers)

write.csv(immune.markers, "immune.markers_vsCD8_Mem_Tcruk.csv")




Tcruk.markers <- FindMarkers(object = reap_REAP, ident.1 = "Tcruk.t",  assay.type = "RNA", logfc.threshold=0.25, only.pos = FALSE)
head(Tcruk.markers)

head(reap_REAP@meta.data)
VlnPlot(reap_REAP, features.plot = c("REAP_CD3","REAP_CD8a","REAP_CD45RA","REAP_CD45RO","REAP_CD27","REAP_CD197", "REAP_CD69","REAP_CD279","REAP_CD25","EOMES","MKI67"), group.by="ident")
quartz.save("violin_tot.pdf",type="PDF")

GenePlot(reap_REAP, gene1="EOMES", gene2="REAP_CD69", cell.ids=WhichCells(reap_REAP,"Tcruk.t"))
quartz.save("GenePlot_CD69_EOMES.pdf",type="PDF")

GenePlot(reap_REAP, gene1="MKI67", gene2="REAP_CD69", cell.ids=WhichCells(reap_REAP,"Tcruk.t"))
quartz.save("GenePlot_CD69_Ki67.pdf",type="PDF")

reap_REAP<-DoKMeans(object= reap_REAP, k.genes=3)
KMeansHeatmap(object= reap_REAP, cex.col=3, cex.row=3, group.cex=3, max.genes=100,assay.type="RNA")
quartz.save("KMeansHeatmap_RNA_Tcruk.pdf",type="PDF")
DoHeatmap(object= reap_REAP, cex.col=3, cex.row=3, group.cex=3, use.scaled=FALSE,assay.type="REAP", col.low="blue", col.high="red")
quartz.save("Heatmap_REAP_Tcruk.pdf",type="PDF")
