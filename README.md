# zebTestis
Codes for the projects of zebrafish testis scRNA-seq analysis. Inlcuding cell ranger preprocess, seurat analysis, GO enrichment analysis, and ligand-receptor analysis.
## fastq2matrix
Convert fastq files into UMI matrix through running cell ranger.
## basicAnalysis 
Analyses including Seurat cell cluster analysis, marker gene finding, and ClusterProfiler GO enrichment analysis.
## LRanalysis
Ligand-receptor analysis for paracrine cell-cell communcation.
## Data availability
[raw sequencing data](https://ngdc.cncb.ac.cn/gsa/browse/CRA003925)

[cellranger matrix files](https://figshare.com/articles/dataset/filtered_feature_bc_matrix_zip/19582615)

[Seurat rds file with cell annotation and UMAP](https://figshare.com/articles/dataset/Seurat_object_with_cell_type_annotation_and_UMAP_coordinates_for_zebrafish_testis_single_cell_RNA_sequencing_datasets/27922725?file=50852895)
## Citation:
Qian P, Kang J, Liu D and Xie G (2022) Single Cell Transcriptome Sequencing of Zebrafish Testis Revealed Novel Spermatogenesis Marker Genes and Stronger Leydig-Germ Cell Paracrine Interactions. Front. Genet. 13:851719. doi: 10.3389/fgene.2022.851719
