library(Seurat)
library(scsrctdb) #SingleCellSignalR based on CellTalkDB
options(future.globals.maxSize = 10000 * 1024^2) # 10G memory
#sample="zebTestis"
#dir="/home/db/private/XieLab/zebrafish_testis_10X/run_cellranger/zeb_testis/outs/filtered_feature_bc_matrix/"
#obj_data <- Read10X(data.dir = dir)
obj=readRDS("../../zebTestis.rds")

#zebrafish-human 1-1 homologs
zeb.hsa=read.table("/home/db/public/Ensembl/EnsemblRelease102/zebrafish/homologs/unique_matched/zebrafish_human_unique.tsv",sep="\t",header=T)
rownames(zeb.hsa)=zeb.hsa$human

obj.matrix=obj@assays$RNA@counts

genes.exp=rownames(obj.matrix)

zeb.gene=zeb.hsa$zebrafish[zeb.hsa$zebrafish %in% genes.exp]
obj.matrix.select=obj.matrix[zeb.gene,]
hsa.gene=zeb.hsa$human[zeb.hsa$zebrafish %in% genes.exp]
#re-name rownames to human homolog gene
rownames(obj.matrix.select)=hsa.gene
clusters=obj$seurat_clusters
clusters.n=as.numeric(as.character(clusters))
obj.matrix.select.s=as.matrix(obj.matrix.select)
#run SingleCellSignalR
cell_signal <- cell_signaling(data = obj.matrix.select.s,
                              genes = hsa.gene,
                              cluster = clusters.n,
			      c.names=unique(clusters.n),
			      s.score=0.1,
                              gene_resive = T,
                              species = 'homo sapiens')

#change back to zebrafish gene name
for(pair_name in names(cell_signal)){
	current_signal=cell_signal[[pair_name]]
	cnames=colnames(current_signal)
	col1_new=zeb.hsa[current_signal[,1],"zebrafish"]
	col2_new=zeb.hsa[current_signal[,2],"zebrafish"]
	nc=ncol(current_signal)
	new_signal=cbind(col1_new,col2_new,current_signal[3:nc])
	colnames(new_signal)=cnames
	cell_signal[[pair_name]]=new_signal
	write.table(new_signal,file=paste0("LR_interactions_",pair_name,"_paracrine_zebGenes.tsv"),sep="\t",quote=F)
}



#visualize_interactions(signal = cell_signal,write.in=c(1,4))
pdf("zebTestis_cell_signaling.pdf")
visualize(cell_signal)
dev.off()
for(i in names(cell_signal)){
	pdf(paste0("zebTestis_cell_signaling_paracrine_show",i,".pdf"))
	visualize(cell_signal,show.in=i)
	dev.off()
}
