#!/bin/bash
fq_dir=../fq/clean/
ref=/home/db/public/SoftwareIndex/cellranger_5.0.0/zebrafish/Danio_rerio
cellranger count --id="zeb_testis" \
                   --transcriptome=$ref \
                   --fastqs=$fq_dir \
                   --sample="ze-testis" \
                   --expect-cells=15000 \
                   --localcores=40 \
                   --localmem=64
