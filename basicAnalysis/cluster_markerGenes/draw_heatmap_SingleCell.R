library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
mycolor=colorRampPalette(c("tomato2","snow1","skyblue2"))(50)
obj = readRDS("../../zebTestis.rds")


cell.data=read.table("../../cell_cluster_ids.tsv",header=F,sep="\t")
colnames(cell.data)=c("cell","cluster")
cells <- cell.data %>% filter(cluster %in% c(3,5)) %>% select(cell)
#standard normalize for RNA assay
#log1p(RPM) stored in RNA assay slot data, natural-log 
#obj <- NormalizeData(obj,verbose = FALSE, normalization.method = "LogNormalize",
#                                          assay="RNA",scale.factor = 1e6) #log1p(RPM)

#obj <- ScaleData(obj,do.scale = TRUE , do.center = TRUE,
#                                  assay="RNA")
obj.markers <- read.table("../SPG_subPopulation_markerGenes.tsv",header=T,sep="\t")
top10 <- obj.markers %>% top_n(n = 10, wt = avg_log2FC)
bottom10 <- obj.markers %>% top_n(n = -10, wt = avg_log2FC)
#top.genes <- obj.markers %>% filter(abs(avg_log2FC) > 1)
genes=as.character(c(top10$gene,bottom10$gene))
#genes=top.genes$gene
pdf("ZebTestis_SPG_TopGenes_Heatmap.pdf",width=5,height=5)
p <- DoHeatmap(obj, features = genes,cells=cells$cell, assay = "RNA",
	       slot = "scale.data", raster=F, disp.min= -2)+ scale_fill_gradientn(colours = rev(mycolor))
#p <- DoHeatmap(obj, features = genes,cells=cells$cell, assay = "RNA",
#	       slot = "data", raster=F)+ scale_fill_gradientn(colours = rev(mycolor))
print(p)
dev.off()
