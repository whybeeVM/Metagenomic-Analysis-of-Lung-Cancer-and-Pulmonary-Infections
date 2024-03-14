#!/bin/bash

set -e
echo "Job Start at `date`"

##############Part3 FeatureCount: bam to expression level###########################

mkdir 04_featureCounts

#The path of the software
featureCounts=~/02.Software/miniconda3/bin/featureCounts

#downloads gtf for the hg38 version of H. spiens (human) from UCSC
hg38_annotation=~/RNAseq/annotation/hg38.refGene.gtf

Hisat2_dir=./03_Hisat2_Mapping
featureCounts_dir=./04_featureCounts
SampleID=SampleID

#featurecount
${featureCounts} -T 4 \
-t exon -g gene_id \
-a ${hg38_annotation} \
-o ${featureCounts_dir}/${SampleID}.hg38_count.txt \
${Hisat2_dir}/${SampleID}.hg38.sort.bam \
1>${featureCounts_dir}/${SampleID}.hg38_count.log 2>&1


#get time end the job
echo "Job finished at:" `date`