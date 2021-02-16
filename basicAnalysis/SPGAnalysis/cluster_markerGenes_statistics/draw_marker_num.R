library(ggplot2)
data=read.table("zebTestis_clusterLevelMarkerGenes_counts.tsv",sep="\t",header=T)
#ClusterID	NumberofMarkerGenes

p=ggplot(data=data,aes(x=Change,y=NumberofMarkerGenes))+
  geom_bar(stat="identity",fill="deepskyblue3")+xlab("")+ylab("# of marker genes")+
  geom_text(aes(label=NumberofMarkerGenes), vjust=1.6, color="white",
	              position = position_dodge(0.9), size=3.5)+

pdf("zebTestis_SPGMarkerGenes_counts.pdf",width=2,height=4)
print(p)
dev.off()
