library("clusterProfiler")
library(ggplot2)
options <- commandArgs(trailingOnly = TRUE)
ont_type=as.character(options[1])
data=read.table("../SPG_subPopulation_markerGenes.tsv",sep="\t",header=T)
genes=as.character(data$gene)
avg_log2FC=data$avg_log2FC
change=sapply(avg_log2FC,function(x){
  if(x>0){
     "UP"
  }else{
     "DOWN"
  }
})
#check keyType
#library(org.Hs.eg.db)
#keytypes(org.Hs.eg.db)
#mouse: org.Mm.eg.db
eg = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Dr.eg.db")
name2id=list()
for(i in c(1:nrow(eg))){
  name=as.character(eg[i,1])
  id=as.character(eg[i,2])
  name2id[[name]]=id
}
gene.f=genes %in% eg$SYMBOL
genes_keep=genes[gene.f]
gene_id=sapply(genes_keep,function(x){name2id[[x]]})
gene_id=as.character(gene_id)
my.df=data.frame(Entrez=gene_id,change=change[gene.f])

fun_type="enrichGO"
formula_res <- compareCluster(Entrez~change, data=my.df, fun=fun_type ,OrgDb="org.Dr.eg.db",ont=ont_type,pvalueCutoff=0.05,pAdjustMethod="BH")
write.table(formula_res,file="GO.tsv",sep="\t")
d=dotplot(formula_res, x=~change, showCategory=6)+ theme(axis.text.x = element_text(angle = 45, hjust=1))
pdf(paste0("zebTestis_GO_",ont_type,".pdf"),width=10,height=6)
print(d)
dev.off()
