library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
mycolor=colorRampPalette(c("skyblue2","tomato2"))(10)
obj = readRDS("../../zebTestis.rds")

#rename
rename.file=read.table("../../clusterRename.tsv",header=T,sep="\t")
names=list()
for(i in c(1:nrow(rename.file))){
        c_id=rename.file[i,1]
        c_name=rename.file[i,2]
        c_name=gsub("_"," ",c_name)
        names[[as.character(c_id)]]=c_name
}

print(names)
clusters=levels(obj)

#rename the clusters
new.names <- sapply(clusters,function(x){names[[as.character(x)]]})
obj <- RenameIdents(obj,new.names)


DefaultAssay(obj)="RNA"
obj.markers <- read.table("../SPG_subPopulation_markerGenes.tsv",header=T,sep="\t")
top_genes <- obj.markers %>% top_n(n = 10, wt = avg_log2FC)
bottom_genes <- obj.markers %>% top_n(n = -10, wt = avg_log2FC)

genes=c(top_genes$gene,bottom_genes$gene)
for(gene in genes){
	try({ 
		#p = FeaturePlot(obj, features = gene, cols = mycolor,slot = "scale.data",ncol=1)
		#pdf(paste0(gene,"_","zebTestis_FeaturePlot.pdf"),height=5,width=5)
		#print(p)
		#dev.off()

		p = VlnPlot(obj,assay = "RNA",idents=c("SPG1","SPG2"), features = gene,slot = "scale.data", ncol=1,pt.size=0)
		pdf(paste0(gene,"_","zebTestis_VlnPlot.pdf"),height=5,width=3.5)
		print(p)
		dev.off()
	})
}
