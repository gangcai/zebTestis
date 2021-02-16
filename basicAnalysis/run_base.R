library(Seurat)
library(ggplot2)
library(sctransform)
suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--min.cell", type="integer", default=3,
    help="min number of cell [default %(default)s]",
    metavar="number")
parser$add_argument("--min.feature", type="integer", default=200,
    help="min number of nFeature_RNA [default %(default)s]",
    metavar="number")
parser$add_argument("--variable.features.n", type="integer", default=3000,
    help="# of variable features [default %(default)s]",
    metavar="number")
parser$add_argument("--max.feature", type="integer", default=5000,
    help="max number of nFeature_RNA [default %(default)s]",
    metavar="number")
parser$add_argument("--pc.num", type="integer", default=30,
    help="number of PCs used [default %(default)s]",
    metavar="number")
parser$add_argument("--regress.type", type="integer", default=1,
    help="types of regress [default %(default)s]",
    metavar="number")
parser$add_argument("--resolution", type="double", default=0.25,
    help="clustering resolution value [default %(default)s]",
    metavar="number")
parser$add_argument("--mt.per", type="integer", default=5,
    help="mt percentage [default %(default)s]",
    metavar="number")
args <- parser$parse_args()

min.cell=args$min.cell
min.feature=args$min.feature
max.feature=args$max.feature
pc.num=args$pc.num
resolution=args$resolution
variable.features.n=args$variable.features.n
regress.type=args$regress.type
mt.per=args$mt.per
parameters=c(min.cell,min.feature,max.feature,pc.num,
             resolution,variable.features.n,regress.type,mt.per)
print(parameters)
suffix=paste(min.cell,min.feature,max.feature,pc.num,
             resolution,variable.features.n,regress.type,mt.per,sep="_")
print(suffix)

options(future.globals.maxSize = 10000 * 1024^2) # 10G memory
sample="zebTestis"
dir="/home/db/private/XieLab/zebrafish_testis_10X/run_cellranger/zeb_testis/outs/filtered_feature_bc_matrix/"
obj_data <- Read10X(data.dir = dir)
obj <- CreateSeuratObject(counts = obj_data,project = sample, min.cells = min.cell, min.features = min.feature)


# store mitochondrial percentage in object meta data
obj <- PercentageFeatureSet(obj, pattern = "^mt-", col.name = "percent.mt")

#QC check before filtering
p0 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
pdf(paste0(sample,"_sctransform_QC_before.pdf"))
print(p0)
dev.off()

obj <- subset(obj, subset = nFeature_RNA > min.feature & nFeature_RNA < max.feature & percent.mt < mt.per)

# run sctransform
if(regress.type == 1){
	print("regress by nFeature_RNA, nCount_RNA, percent.mt")
	obj <- SCTransform(obj,variable.features.n=variable.features.n,
			   vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"), verbose = FALSE)
}
if(regress.type == 2){
	print("regress by nFeature_RNA, nCount_RNA")
	obj <- SCTransform(obj,variable.features.n=variable.features.n,
			   vars.to.regress = c("nFeature_RNA", "nCount_RNA"), verbose = FALSE)
}
if(regress.type == 3){
	print("regress by nCount_RNA")
	obj <- SCTransform(obj,variable.features.n=variable.features.n,
			   vars.to.regress = c("nCount_RNA"), verbose = FALSE)
}


# These are now standard steps in the Seurat workflow for visualization and clustering
obj <- RunPCA(obj, verbose = FALSE, npcs = pc.num)
obj <- RunUMAP(obj, dims = 1:pc.num, verbose = FALSE)
obj <- RunTSNE(obj, dims = 1:pc.num, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:pc.num, verbose = FALSE)
obj <- FindClusters(obj,resolution = resolution,  verbose = FALSE)

p <- DimPlot(obj, label = TRUE, reduction="tsne") + NoLegend()
pdf(paste0(suffix,"_",sample,"_Cluster_tSNE_Plot.pdf"))
print(p)
dev.off()

p <- DimPlot(obj, label = TRUE, reduction="umap") + NoLegend()
pdf(paste0(suffix,"_",sample,"_Cluster_UMAP_Plot.pdf"))
print(p)
dev.off()

#standard normalize for RNA assay
#log1p(RPM) stored in RNA assay slot data, natural-log 
obj <- NormalizeData(obj,verbose = FALSE, normalization.method = "LogNormalize",
                                          assay="RNA",scale.factor = 1e6) #log1p(RPM)

obj <- ScaleData(obj,do.scale = TRUE , do.center = TRUE,
                                  assay="RNA")

saveRDS(obj, file = paste0(sample,".rds"))

# find markers for every cluster compared to all remaining cells, report only the positive ones
obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(obj.markers,file=paste0(sample,"_marker_genes.tsv"),sep="\t",quote=F,row.names=F)


#QC plot
p2 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
pdf("testis_sctransform_QC_After.pdf")
print(p2)
dev.off()

clusters=obj$seurat_clusters
clusters=data.frame(clusters)
write.table(clusters,file="cell_cluster_ids.tsv",sep="\t",row.names=T,col.names=F,quote=F)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

