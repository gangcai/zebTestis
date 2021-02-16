import re
outfile=open("zebTestis_clusterLevelMarkerGenes_counts.tsv","w")
outfile.write("Change\tNumberofMarkerGenes\n")
i=0
counts={"UP":0,"DOWN":0}
for line in open("../SPG_subPopulation_markerGenes.tsv"):
	i+=1
	if i==1:
		continue
	items=line.rstrip().split("\t")
	cid="c"+items[-2]
	fc=float(items[2])
	if fc > 0:
		counts["UP"]+=1
	else:
		counts["DOWN"]+=1




for cid in counts.keys():
	n=counts[cid]
	outfile.write(cid+"\t"+str(n)+"\n")
outfile.close()
