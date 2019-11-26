library(readxl)
library(matrixStats)
library(flowCore)
library(FlowSOM)
library(ComplexHeatmap)
library(cowplot)
library(lme4)
library(multcomp)



rm(list=ls())

setwd("/.../FlowRepository_FR-FCM-ZYCC_files")
metadata_filename<-"metadata_M.xlsx"

md<-read_excel(metadata_filename)
md$condition<-factor(md$condition, levels=c("M"))
color_conditions <- c("#009E73")#"#6A3D9A", "#FF7F00","#CC79A7", 
names(color_conditions) <- levels(md$condition)


#fcs_filename <- "FlowRepository_FR-FCM-ZY34_files.zip.download"
#setwd("/Users/svalpione/Downloads/FlowRepository_FR-FCM-ZY34_files")
#download.file(paste0("/Users/svalpione/Downloads", "/", fcs_filename), destfile = fcs_filename, mode = "wb") 
#unzip(fcs_filename) 
#fcs_raw <- read.flowSet(md$file_name, transformation = FALSE, truncate_max_range = FALSE)

panel_filename<-"MC_panel3.xlsx"
panel <- read_excel(panel_filename)

head(data.frame(md))
library(flowCore)

fcs_raw <- read.flowSet(md$filename, transformation = FALSE, truncate_max_range = FALSE, invert.pattern=TRUE, column.pattern="SampleID", alter.names=TRUE)
fcs_raw[[1]]
traceback()




head(data.frame(panel)) 
panel$Antigen <- gsub("-", ".", panel$Antigen)
#panel$Antigen[panel$Antigen=="CD4_CD4"]<-"CD4"
panel$Antigen1<-paste(panel$Antigen, panel$sep, sep="")
panel$Antigen2<-paste(panel$Antigen1, panel$Isotope, sep="") 
panel$Antigen2 <- gsub("/", "", panel$Antigen2)
panel$Antigen<-panel$Antigen2

panel_fcs <- pData(parameters(fcs_raw[[1]]))
head(panel_fcs$desc)
panel_fcs<-panel_fcs[!(panel_fcs$name=="SampleID"),]#check if you need this
class(panel_fcs)



panel_fcs$desc<-panel_fcs$name
#panel_fcs$desc <- gsub("-", "_", panel_fcs$desc)
#panel_fcs$desc<-paste("X", panel_fcs$desc,sep="")
#panel_fcs$desc[panel_fcs$desc=="X145Nd_CD4.CD4"]<-"X145Nd_CD4"
lineage_markers <- panel$Antigen[panel$T_cell == 1]
functional_markers <- panel$Antigen[panel$Functional == 1]
Tcell<-panel$Antigen[panel$T_cell==1]

class(panel_fcs$desc)

all(lineage_markers %in% panel_fcs$desc)
all(functional_markers %in% panel_fcs$desc)
which(functional_markers %in% panel_fcs$desc)


which(Tcell %in% panel_fcs$desc)

all(Tcell %in% panel_fcs$desc)#CD8 in 177 or 168?!?
#"170Er_CD3" %in% panel_fcs$desc #his is why it doesn't work!

head(fcs_raw)

fcs <- fsApply(fcs_raw, function(x, cofactor=5){  #try with cofactor=150
  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  expr <- asinh(expr [, c(Tcell)]/ cofactor) #try with c(lineage_markers, functional_markers)  
  exprs(x) <- expr
  x
})
fcs 

expr <- fsApply(fcs, exprs)
dim(expr)
library(matrixStats)
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1




sample_ids <- rep(md$sample_id, fsApply(fcs, nrow))
batch<-rep(md$batch, fsApply(fcs, nrow))
batch<-as.matrix(batch)
colnames(batch)<-paste("batch")
condition<-rep(md$condition, fsApply(fcs, nrow))


#library(limma)
#Texpr01<-t(expr01)
#new<-t(removeBatchEffect(Texpr01,batch=batch))
#expr01<-new


library(rms)
describe(expr01)
#colnames(new)<-colnames(expr01)

head(expr01)

sample_ids <- rep(md$sample_id, fsApply(fcs, nrow))

library(ggplot2)
library(reshape2)

ggdf <- data.frame(sample_id = sample_ids, expr)
ggdf <- melt(ggdf, id.var = "sample_id",
  value.name = "expression", variable.name = "antigen")
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = expression, color = condition, group = sample_id)) + geom_density() + facet_wrap(~ antigen, nrow = 4, scales = "free") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.text = element_text(size = 7), axis.text = element_text(size = 5)) + scale_color_manual(values = color_conditions)

library(dplyr)

cell_table <- table(sample_ids)
ggdf <- data.frame(sample_id = names(cell_table), cell_counts = as.numeric(cell_table))
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = condition)) + geom_bar(stat = "identity")+ geom_text(aes(label = cell_counts), hjust=0.5, vjust=-0.5, size = 2.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_fill_manual(values = color_conditions, drop = FALSE) + scale_x_discrete(drop = FALSE) 


expr_median_sample_tbl <- data.frame(sample_id = sample_ids, expr) %>% group_by(sample_id) %>%  summarize_all(funs(median))
expr_median_sample <- t(expr_median_sample_tbl[, -1])
colnames(expr_median_sample) <- expr_median_sample_tbl$sample_id

library(limma)

mds <- plotMDS(expr_median_sample, plot = FALSE)

library(ggrepel)

ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
  sample_id = colnames(expr_median_sample))
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = MDS1, y = MDS2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = sample_id)) +
  theme_bw() +
scale_color_manual(values = color_conditions) + coord_fixed() 
quartz.save("batch_effect_Melanomas.pdf", type="pdf")


library(RColorBrewer)
library(pheatmap)

mm <- match(colnames(expr_median_sample), md$sample_id)
annotation_col <- data.frame(condition = md$condition[mm], row.names = colnames(expr_median_sample))
annotation_colors <- list(condition = color_conditions)

color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)

pheatmap(expr_median_sample, color = color, display_numbers = TRUE, number_color = "black", fontsize_number = 5, annotation_col = annotation_col, annotation_colors = annotation_colors, clustering_method = "average") 





NRS <- function(x, ncomp = 3){
	pr <- prcomp(x, center = TRUE, scale. = FALSE)
	score <- rowSums(outer(rep(1, ncol(x)),
	pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
	return(score)
	}

	
nrs_sample <- fsApply(fcs[, Tcell], NRS, use.exprs = TRUE)

rownames(nrs_sample) <- md$sample_id
nrs <- colMeans(nrs_sample, na.rm = TRUE)

lineage_markers_ord <- names(sort(nrs, decreasing = TRUE))
nrs_sample <- data.frame(nrs_sample, check.names = FALSE)#set check.names=FALSE to avoid X
nrs_sample$sample_id <- rownames(nrs_sample)

ggdf <- melt(nrs_sample, id.var = "sample_id", value.name = "nrs", variable.name = "antigen")

head(nrs_sample)

ggdf$antigen <- factor(ggdf$antigen, levels = lineage_markers_ord)
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = antigen, y = nrs)) + geom_point(aes(color = condition), alpha = 0.9, position = position_jitter(width = 0.3, height = 0)) + geom_boxplot(outlier.color = NA, fill = NA)+ stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_conditions) 

library(FlowSOM) 
fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE) 
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_markers)
#Get the cell clustering into 100 SOM codes
cell_clustering_som <- som$map$mapping[,1]

library(ConsensusClusterPlus)
codes <- som$map$codes
plot_outdir <- "consensus_plots"
nmc <- 20

mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100, pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png", clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average", distance = "euclidean", seed = 1234)
 
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[cell_clustering_som]

color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
#head(mc)


names<-c("CD45RO","CD45RA","CD8","CD3","CD4","CD27","CCR7")#the order is checked using the heatmap using expr01, be VERY carefull to double check the final order of the heatmap columns before changing colnames(expr_heat)<-names

plot_clustering_heatmap_wrapper <- function(expr, expr01, cell_clustering, color_clusters, cluster_merging = NULL){
	expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>% group_by(cell_clustering) %>% summarize_all(funs(median))
	expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>% group_by(cell_clustering) %>% summarize_all(funs(median))
	clustering_table <- as.numeric(table(cell_clustering))
	clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
	d <- dist(expr_median[, colnames(expr)], method = "euclidean")
	cluster_rows <- hclust(d, method = "average")
	expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
	colnames(expr_heat) <- names #added
	rownames(expr_heat) <- expr01_median$cell_clustering
	color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100) 
	legend_breaks = seq(from = 0, to = 1, by = 0.2)
	labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop ,"%)") 
	annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering))
	rownames(annotation_row) <- rownames(expr_heat)
	color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)] 
	names(color_clusters1) <- levels(annotation_row$Cluster) 
	annotation_colors <- list(Cluster = color_clusters1)
	if(!is.null(cluster_merging)){
	cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
	annotation_row$Cluster_merging <- cluster_merging$new_cluster
	color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
	names(color_clusters2) <- levels(cluster_merging$new_cluster)
	annotation_colors$Cluster_merging <- color_clusters2
	}
	pheatmap(expr_heat, color = color_heat, cluster_cols = FALSE,
	cluster_rows = cluster_rows, labels_row = labels_row,
	display_numbers = TRUE, number_color = "black",
	fontsize = 6, fontsize_number = 6,  legend_breaks = legend_breaks,
	annotation_row = annotation_row, annotation_colors = annotation_colors)
	}

#cluster 9 is Tcruks in tot
I<-plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord], expr01 = expr01[,lineage_markers_ord], cell_clustering = cell_clustering1, color_clusters = color_clusters)
I
quartz.save("plot_clustering_heatmap_wrapper.pdf", type="pdf")

#median expression data

library(ggridges)

plot_clustering_distr_wrapper <- function(expr, cell_clustering){
cell_clustering <- factor(cell_clustering)
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
  group_by(cell_clustering) %>% summarize_all(funs(median))
d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
freq_clust <- table(cell_clustering)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  cell_clustering <- factor(cell_clustering,
    labels = paste0(levels(cell_clustering), " (", freq_clust, "%)"))
  
ggd <- melt(data.frame(cluster = cell_clustering, expr), id.vars = "cluster", value.name = "expression",variable.name = "antigen")
ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
ggd$reference <- "no"
### The reference data
ggd_bg <- ggd
ggd_bg$cluster <- "reference"
ggd_bg$reference <- "yes"
ggd_plot <- rbind(ggd, ggd_bg)
ggd_plot$cluster <- factor(ggd_plot$cluster, levels = c(levels(cell_clustering)[rev(cluster_rows$order)], "reference"))

ggplot() + geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster, color = reference, fill = reference), alpha = 0.3) + facet_wrap( ~ antigen, scales = "free_x", nrow = 2) +    
theme_ridges() + theme(axis.text = element_text(size = 7), strip.text = element_text(size = 7), legend.position = "none")
}

plot_clustering_distr_wrapper(expr = expr[, Tcell], cell_clustering = cell_clustering1)
quartz.save("cluster_median_expression.pdf", type="pdf")

library(ComplexHeatmap)

#t-SNE
  
dups <- which(!duplicated(expr[, Tcell]))
inds <- split(1:length(sample_ids), sample_ids)
tsne_ncells <- pmin(table(sample_ids), 2000)
set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})

tsne_inds <- unlist(tsne_inds)
tsne_expr <- expr[tsne_inds, Tcell]

library(Rtsne)
set.seed(1234)

#occhio che questo richiede parecchi minuti
tsne_out <- Rtsne(tsne_expr, check_duplicates = FALSE, pca = FALSE)

dr <- data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2], expr[tsne_inds, Tcell])

par(mfrow=c(2,4))
A<-ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = CD3.170)) + geom_point(size = 0.8) + theme_bw() + scale_color_gradientn("CD3", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)) 
B<-ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = CD8.177)) + geom_point(size = 0.8) + theme_bw() + scale_color_gradientn("CD8", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)) 
C<-ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = CD4.145)) + geom_point(size = 0.8) + theme_bw() + scale_color_gradientn("CD4", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)) 
D<-ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = CD45RA.153)) + geom_point(size = 0.8) + theme_bw() + scale_color_gradientn("CD45RA", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)) 
E<-ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = CD45RO.164)) + geom_point(size = 0.8) + theme_bw() + scale_color_gradientn("CD45RO", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)) 
F<-ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = CD27.167)) + geom_point(size = 0.8) + theme_bw() + scale_color_gradientn("CD27", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)) 
G<-ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = CCR7.159)) + geom_point(size = 0.8) + theme_bw() + scale_color_gradientn("CCR7", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50)) 


#quartz.save("t-SNE.pdf", type="pdf")


dr$sample_id <- sample_ids[tsne_inds]
mm <- match(dr$sample_id, md$sample_id)
dr$condition <- md$condition[mm]
dr$cell_clustering1 <- factor(cell_clustering1[tsne_inds], levels = 1:nmc)

H <- ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering1)) + geom_point(size = 0.8) +theme_bw() + scale_color_manual(values = color_clusters) + guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
H
quartz.save("t-SNE_clusters.pdf", type="pdf")




#clusters with Tcruk name for cluster 9

cluster_merging1 <- data.frame("original_cluster"=c(1:20), "new_cluster"=c(1:20)) 
cluster_merging1$new_cluster[cluster_merging1$new_cluster=="9"]<-"Tcruk"
data.frame(cluster_merging1)
levels_clusters_merged <- c("1","2","3","4","5","6","7","8","Tcruk","10","11","12","13","14","15","16","17","18","19","20")
cluster_merging1$new_cluster <- factor(cluster_merging1$new_cluster, levels = levels_clusters_merged)

mm <- match(cell_clustering1, cluster_merging1$original_cluster) 
cell_clustering1m <- cluster_merging1$new_cluster[mm]
mm <- match(code_clustering1, cluster_merging1$original_cluster)
code_clustering1m <- cluster_merging1$new_cluster[mm]
dr$cell_clustering1m <- cell_clustering1m[tsne_inds]

H<-ggplot(dr,  aes(x = tSNE1, y = tSNE2, color = cell_clustering1m)) + geom_point(size = 0.8) + theme_bw() + scale_color_manual(values = color_clusters) + guides(color = guide_legend(override.aes = list(size = 4), ncol=2))
quartz.save("t-SNE_Tcruk.pdf", type="pdf")



plot_clustering_heatmap_wrapper <- function(expr, expr01, cell_clustering, color_clusters, cluster_merging = NULL){
	expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>% group_by(cell_clustering) %>% summarize_all(funs(median))
	expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>% group_by(cell_clustering) %>% summarize_all(funs(median))
	clustering_table <- as.numeric(table(cell_clustering))
	clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
	d <- dist(expr_median[, colnames(expr)], method = "euclidean")
	cluster_rows <- hclust(d, method = "average")
	expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
	colnames(expr_heat) <- names #added
	expr01_median$cluster_merging1<-c("1","2","3","4","5","6","7","8","Tcruk","10","11","12","13","14","15","16","17","18","19","20")#added
	rownames(expr_heat) <- expr01_median$cluster_merging1 # it was cell_clustering
	color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100) 
	legend_breaks = seq(from = 0, to = 1, by = 0.2)
	labels_row <- paste0(expr01_median$cluster_merging1, " (", clustering_prop ,"%)") 
	annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering))
	rownames(annotation_row) <- rownames(expr_heat)
	color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)] 
	names(color_clusters1) <- levels(annotation_row$Cluster) 
	annotation_colors <- list(Cluster = color_clusters1)
	if(!is.null(cluster_merging)){
	cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
	annotation_row$Cluster_merging <- cluster_merging$new_cluster
	color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
	names(color_clusters2) <- levels(cluster_merging$new_cluster)
	annotation_colors$Cluster_merging <- color_clusters2
	}
	pheatmap(expr_heat, color = color_heat, cluster_cols = FALSE,
	cluster_rows = cluster_rows, labels_row = labels_row,
	display_numbers = TRUE, number_color = "black",
	fontsize = 6, fontsize_number = 6,  legend_breaks = legend_breaks,
	annotation_row = annotation_row, annotation_colors = annotation_colors)
	}

#cluster 9 is Tcruks in tot
I<-plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord], expr01 = expr01[,lineage_markers_ord], cell_clustering = cell_clustering1, color_clusters = color_clusters)
I



library(grid)
library(ggplotify)
I<-as.grob(I)

library(cowplot)
plot_grid(A, B,C,D,E,F,G,H, I, labels="AUTO")
quartz.save("t-SNE_tumor_melanoma.pdf", type="pdf")
dev.off()








library(lme4)
library(multcomp)


model.matrix( ~ condition_time, data = md)#model formula without random effect
contrast_names <- c("all")
k1 <- c(0, 1)
K <- matrix(k1, nrow = 1, byrow = TRUE, dimnames = list(contrast_names))
K

FDR_cutoff <- 0.05
counts_table <- table(cell_clustering1m, sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

ggdf <- melt(data.frame(cluster = rownames(props), props), id.vars = "cluster", value.name = "proportion", variable.name = "sample_id") 

ggdf$cluster <- factor(ggdf$cluster, levels = levels_clusters_merged) ## Add condition info

ggdf$sample_id<-gsub("X", "", ggdf$sample_id)


mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- factor(md$condition_time[mm])

library(forcats)
color_conditions1 <- c("#009E73","#6A3D9A")


fct_explicit_na (ggdf$condition, "NA") 
names(color_conditions1) <- levels(md$condition_time)


head(md)

ggplot(ggdf, aes(x = sample_id, y = proportion, fill = cluster)) + geom_bar(stat = "identity") + facet_wrap(~condition, scales = "free_x") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = color_clusters)
quartz.save("clusters in samples.pdf", type="PDF")
ggdf$patient_id <- factor(md$patient_id[mm])

md$condition[mm]

library(forcats)
color_conditions1 <- c("#009E73","#6A3D9A", "#FF7F00","#CC79A7", "#A89DA3")
fct_explicit_na(md$condition_tot, "NA")
names(color_conditions1) <- levels(md$condition_tot)

ggplot(ggdf) + geom_boxplot(aes(x = condition, y = proportion, color = condition, fill = condition), position = position_dodge(), alpha = 0.5, outlier.color = NA) + geom_point(aes(x = condition, y = proportion, color = condition, shape = patient_id), alpha = 0.8, position = position_jitterdodge()) + facet_wrap(~ cluster, scales = "free", nrow = 2) + theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),  axis.title.x = element_blank(), strip.text = element_text(size = 6)) + scale_color_manual(values = color_conditions1) + scale_fill_manual(values = color_conditions1) + scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) 
quartz.save("clusters over treatment_BRAFi.pdf", type="PDF")

dev.off()



library(multcomp)
model.matrix( ~ condition_time, data = md)


contrast_names <- c("BasevsTx")
k1 <- c(0,1)
K <- matrix(k1, nrow = 1, byrow = TRUE, dimnames = list(contrast_names))
K

FDR_cutoff <- 0.05
counts_table <- table(cell_clustering1m, sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

md$sample_id<-as.character(md$sample_id)#this is necessary if sample_id is numeric in the metadata

formula_glmer_binomial1 <- y/total ~ condition_time+ (1|sample_id) 
formula_glmer_binomial2 <- y/total ~ condition_time + (1|patient_id) + (1|sample_id) 

md$condition_time

differential_abundance_wrapper <- function(counts, md, formula, K){ 
	ntot <- colSums(counts)
	fit_binomial <- lapply(1:nrow(counts), function(i){ 
		data_tmp <- data.frame(y = as.numeric(counts[i, md$sample_id]), total = ntot[md$sample_id], md)
		fit_tmp <- glmer(formula, weights = total, family = binomial, data = data_tmp)
		 
		out <- apply(K, 1, function(k){
			contr_tmp <- glht(fit_tmp, linfct = matrix(k,1)) 
			summ_tmp <- summary(contr_tmp)
			pval <- summ_tmp$test$pvalues
			return(pval)
		})
		return(out)
	})
	pvals <- do.call(rbind, fit_binomial) 
	colnames(pvals) <- paste0("pval_", contrast_names)
	rownames(pvals) <- rownames(counts)
	## Adjust the p-values
	adjp <- apply(pvals, 2, p.adjust, method = "BH") 
	colnames(adjp) <- paste0("adjp_", contrast_names) 
	return(list(pvals = pvals, adjp = adjp))
}

da_out1 <- differential_abundance_wrapper(counts, md = md, formula = formula_glmer_binomial1, K = K)
apply(da_out1$adjp < FDR_cutoff, 2, table)



da_out2 <- differential_abundance_wrapper(counts, md = md, formula = formula_glmer_binomial2, K = K)
apply(da_out2$adjp < FDR_cutoff, 2, table)
print(da_out2)


da_output2 <- data.frame(cluster = rownames(props), props, da_out2$pvals, da_out2$adjp, row.names = NULL)
print(head(da_output2), digits = 2)
write.csv(da_output2, "da_output2_tot_cof5.csv")

















############

#END of analyses as per 9Apr19

######   NO - begin


normalization_wrapper <- function(expr, th = 2.5){
	expr_norm <- apply(expr, 1, function(x){
		sdx <- sd(x, na.rm = TRUE)
		if(sdx == 0){
			x <- (x - mean(x, na.rm = TRUE))
		}else{
			x <- (x - mean(x, na.rm = TRUE)) / sdx
		}
		x[x > th] <- th 
		x[x < -th] <- -th
		return(x)
	})
	expr_norm<-t(expr_norm)
}

plot_differential_heatmap_wrapper <- function(expr_norm, sign_adjp, condition_time, color_conditions, th = 2.5){
	oo <- order(condition_time)
	condition_time <- condition_time[oo]
	expr_norm <- expr_norm[, oo, drop = FALSE]
	labels_row <- paste0(rownames(expr_norm), " (", sprintf( "%.02e", sign_adjp), ")")
	labels_col <- colnames(expr_norm)
	annotation_col <- data.frame(condition_time = factor(condition_time))
	rownames(annotation_col) <- colnames(expr_norm)
	annotation_colors <- list(condition = color_conditions)
	color <- colorRampPalette(c("#87CEFA", "#56B4E9", "#56B4E9", "#0072B2",
	"#000000", "#D55E00", "#E69F00", "#E69F00", "#FFD700"))(100)
	breaks = seq(from = -th, to = th, length.out = 101)
	legend_breaks = seq(from = -round(th), to = round(th), by = 1)
	gaps_col <- as.numeric(table(annotation_col$condition_time))
	pheatmap(expr_norm, color = color, breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, annotation_col = annotation_col, annotation_colors = annotation_colors, fontsize = 2)
} 






asin_table <- asin(sqrt((t(t(counts_table) / colSums(counts_table)))))
asin <- as.data.frame.matrix(asin_table)
sign_clusters <- names(which(sort(da_out2$adjp[, "adjp_BasevsTx"]) < FDR_cutoff)) 
sign_adjp <- da_out2$adjp[sign_clusters , "adjp_BasevsTx", drop=FALSE]
## Normalize the transformed proportions to mean = 0 and sd = 1
asin_norm <- normalization_wrapper(asin[sign_clusters, ])

mm <- match(colnames(asin_norm), md$sample_id) 

plot_differential_heatmap_wrapper(expr_norm = asin_norm, sign_adjp = sign_adjp, condition = md$condition_time[mm], color_conditions = color_conditions)


###### NO - END




fcs <- fsApply(fcs_raw, function(x, cofactor=5){  #try with cofactor=150
  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  expr <- asinh(expr [, c(Tcell, functional_markers)]/ cofactor) #try with c(lineage_markers, functional_markers)  
  exprs(x) <- expr
  x
})
fcs 






expr <- fsApply(fcs, exprs)
dim(expr)
library(matrixStats)
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

expr_median_sample_cluster_tbl <- data.frame(expr[, functional_markers], sample_id = sample_ids, cluster = cell_clustering1m) %>% group_by(sample_id, cluster) %>% summarize_all(funs(median))

expr_median_sample_cluster_melt <- melt(expr_median_sample_cluster_tbl, id.vars = c("sample_id", "cluster"), value.name = "median_expression", variable.name = "antigen")

expr_median_sample_cluster <- dcast(expr_median_sample_cluster_melt, cluster + antigen ~ sample_id, value.var = "median_expression") 
rownames(expr_median_sample_cluster) <- paste0(expr_median_sample_cluster$cluster, "_", expr_median_sample_cluster$antigen) 
clusters_keep <- names(which((rowSums(counts < 5) == 0)))
keepLF <- expr_median_sample_cluster$cluster %in% clusters_keep 
expr_median_sample_cluster <- expr_median_sample_cluster[keepLF, ]
keep0 <- rowSums(expr_median_sample_cluster[, md$sample_id]) > 0
expr_median_sample_cluster <- expr_median_sample_cluster[keep0, ]

ggdf <- expr_median_sample_cluster_melt[expr_median_sample_cluster_melt$cluster %in% clusters_keep, ] 
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- factor(md$condition[mm])
ggdf$patient_id <- factor(md$patient_id[mm])
ggplot(ggdf) + geom_boxplot(aes(x = antigen, y = median_expression, color = condition, fill = condition), position = position_dodge(), alpha = 0.5, outlier.color = NA) +  geom_point(aes(x = antigen, y = median_expression, color = condition, shape = patient_id), alpha = 0.8, position = position_jitterdodge(), size = 0.7) + facet_wrap(~ cluster, scales = "free_y", ncol=2) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_color_manual(values = color_conditions) + scale_fill_manual(values = color_conditions) + scale_shape_manual(values = c(1:20)) + guides(shape = guide_legend(override.aes = list(size = 2))) 
quartz.save("differential_expression_single_markers_8Jan.pdf", type="pdf")


differential_expression_wrapper <- function(expr_median, md, model = "lmer", formula, K){
	fit_gaussian <- lapply(1:nrow(expr_median), function(i){ 
	data_tmp <- data.frame(y = as.numeric(expr_median[i, md$sample_id]), md)
	switch(model,
	lmer = {
		fit_tmp <- lmer(formula, data = data_tmp)
		}, 
		lm = {
			fit_tmp <- lm(formula, data = data_tmp)
		}) 
	out <- apply(K, 1, function(k){
		contr_tmp <- glht(fit_tmp, linfct = matrix(k,1)) 
		summ_tmp <- summary(contr_tmp)
		pval <- summ_tmp$test$pvalues 
		return(pval)
	})
	return(out)
})
pvals <- do.call(rbind, fit_gaussian) 
colnames(pvals) <- paste0("pval_", contrast_names) 
rownames(pvals) <- rownames(expr_median)
adjp <- apply(pvals, 2, p.adjust, method = "BH") 
colnames(adjp) <- paste0("adjp_", contrast_names) 
return(list(pvals = pvals, adjp = adjp))
}

formula_lm <- y ~ condition_resp +condition_time
formula_lmer <- y ~ condition_resp +condition_time + (1|patient_id) 
de_out1 <- differential_expression_wrapper(expr_median = expr_median_sample_cluster, md = md, model = "lm", formula = formula_lm, K = K) 
apply(de_out1$adjp < FDR_cutoff, 2, table)

de_out2 <- differential_expression_wrapper(expr_median = expr_median_sample_cluster, md = md, model = "lmer", formula = formula_lmer, K = K) 
apply(de_out2$adjp < FDR_cutoff, 2, table)

de_output2 <- data.frame(expr_median_sample_cluster, de_out2$pvals, de_out2$adjp, row.names = NULL)
print(head(de_output2), digits = 2)
write.csv(de_output2,"de_output2func_tot_Jan18.csv")




sign_clusters_markers <- names(which(de_out2$adjp[, "adjp_tx-NR"] < FDR_cutoff)) 
oo <- order(expr_median_sample_cluster[sign_clusters_markers, "cluster"], de_out2$adjp[sign_clusters_markers, "adjp_tx-NR"]) 
sign_clusters_markers <- sign_clusters_markers[oo] 
sign_adjp <- de_out2$adjp[sign_clusters_markers , "adjp_base-NR"] 
expr_s <- expr_median_sample_cluster[sign_clusters_markers,md$sample_id]
expr_median_sample_cluster_norm <- normalization_wrapper(expr_s)
mm <- match(colnames(expr_median_sample_cluster_norm), md$sample_id)



plot_differential_heatmap_wrapper <- function(expr_norm, sign_adjp, condition, color_conditions, th = 2.5){
	oo <- order(condition)
	condition <- condition[oo]
	expr_norm <- expr_norm[, oo, drop = FALSE]
	labels_row <- paste0(rownames(expr_norm), " (", sprintf( "%.02e", sign_adjp), ")")
	labels_col <- colnames(expr_norm)
	annotation_col <- data.frame(condition = factor(condition))
	rownames(annotation_col) <- colnames(expr_norm)
	annotation_colors <- list(condition = color_conditions)
	color <- colorRampPalette(c("#87CEFA", "#56B4E9", "#56B4E9", "#0072B2",
	"#000000", "#D55E00", "#E69F00", "#E69F00", "#FFD700"))(100)
	breaks = seq(from = -th, to = th, length.out = 101)
	legend_breaks = seq(from = -round(th), to = round(th), by = 1)
	gaps_col <- as.numeric(table(annotation_col$condition))
	pheatmap(expr_norm, color = color, breaks = breaks, legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, annotation_col = annotation_col, annotation_colors = annotation_colors, fontsize = 2.5)
} 


plot_differential_heatmap_wrapper(expr_norm = expr_median_sample_cluster_norm, sign_adjp = sign_adjp, condition = md$condition[mm], color_conditions = color_conditions)

quartz.save("differential_heatmap_wrapper_tot.pdf", type="pdf")
#not working, I think you need to go back and include functional_markers in ggdf
ggdf <- melt(data.frame(expr_median_sample[functional_markers, ], antigen = functional_markers), id.vars = "antigen", value.name = "median_expression", variable.name = "sample_id")







 














 






traceback()

colnames(expr)


head(color_clusters)


library(ggridges)
plot_clustering_distr_wrapper <- function(expr, cell_clustering){
	cell_clustering <- factor(cell_clustering)
	expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>% group_by(cell_clustering) %>% summarize_all(funs(median))
	d <- dist(expr_median[, colnames(expr)], method = "euclidean")
	cluster_rows <- hclust(d, method = "average")
	freq_clust <- table(cell_clustering)
	freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
	cell_clustering <- factor(cell_clustering,
	labels = paste0(levels(cell_clustering), " (", freq_clust, "%)"))
	ggd <- melt(data.frame(cluster = cell_clustering, expr), id.vars = "cluster", value.name = "expression", variable.name = "antigen")
	ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
	ggd$reference <- "no"
	ggd_bg <- ggd
	ggd_bg$cluster <- "reference"
	ggd_bg$reference <- "yes"
	ggd_plot <- rbind(ggd, ggd_bg)
	ggd_plot$cluster <- factor(ggd_plot$cluster, levels = c(levels(cell_clustering [rev(cluster_rows$order)], "reference"))

ggplot() + geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster, color = reference, fill = reference), alpha = 0.3) + facet_wrap( ~ antigen, scales = "free_x", nrow = 2) + theme_ridges() + theme(axis.text = element_text(size = 7), strip.text = element_text(size = 7), legend.position = "none")

}

plot_clustering_distr_wrapper(expr = expr[, lineage_markers_ord], cell_clustering = cell_clustering1)



########################
library(ComplexHeatmap)
plot_clustering_heatmap_wrapper2 <- function(expr, expr01,
  Tcell, functional_markers = NULL, sample_ids = NULL,
  cell_clustering, color_clusters, cluster_merging = NULL,
  plot_cluster_annotation = TRUE){
  	
  	expr_median <- data.frame(expr[, Tcell],
    cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  expr01_median <- data.frame(expr01[, Tcell],
    cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
    
clustering_table <- as.numeric(table(cell_clustering))
clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2) 
d <- dist(expr_median[, Tcell], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, lineage_markers])


expr_median_sample_cluster_tbl <- data.frame(expr01[, functional_markers,
    drop = FALSE], sample_id = sample_ids, cluster = cell_clustering) %>%
    group_by(sample_id, cluster) %>% summarize_all(funs(median))

color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100) 
legend_breaks = seq(from = 0, to = 1, by = 0.2)
labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop , "%)") 
annotation_row1 <- data.frame(Cluster = factor(expr01_median$cell_clustering))
color_clusters1 <- color_clusters[1:nlevels(annotation_row1$Cluster)] 
names(color_clusters1) <- levels(annotation_row1$Cluster) 
if(!is.null(cluster_merging)){
	mm <- match(annotation_row1$Cluster, cluster_merging$original_cluster)
	annotation_row2 <- data.frame(Cluster_merging = factor(cluster_merging$new_cluster[mm]))
	color_clusters2 <- color_clusters[1:nlevels(annotation_row2$Cluster_merging)] 
  }


ha1 <- Heatmap(annotation_row1, name = "Cluster",
col = color_clusters1, cluster_columns = FALSE, 
cluster_rows = cluster_rows, row_dend_reorder = FALSE, 
show_row_names = FALSE, 
width = unit(0.5, "cm"), 
rect_gp = gpar(col = "grey")) 
if(!is.null(cluster_merging)){
	ha2 <- Heatmap(annotation_row2, 
name = "Cluster \nmerging", 
col = color_clusters2, 
cluster_columns = FALSE, 
cluster_rows = cluster_rows, 
row_dend_reorder = FALSE, 
show_row_names = FALSE, 
width = unit(0.5, "cm")
}






