#!/bin/bash

set -e
echo "Job Start at `date`"

##############Part2 Mapping: cleandata map to hg38###########################

mkdir 03_Hisat2_Mapping

#The path of the software
hisat2=~/02.Software/hisat2-2.2.1/hisat2
samtools=~/02.Software/miniconda3/bin/samtools

#downloads hisat2 index for the hg38 version of H. spiens (human) from hisat2 website
hisat2_hg38_index=~/RNAseq/hisat2_index/hg38/genome

CleanData_dir=./02_CleanData
Hisat2_dir=./03_Hisat2_Mapping
SampleID=SampleID

#hisat2 mapping
${hisat2} -p 8 --dta \
-x ${hisat2_hg38_index} \
-U ${CleanData_dir}/${SampleID}.clean.fastq.gz  | \
${samtools} view --threads 8 -bS | \
${samtools} sort --threads 8 -m 4G \
-o ${Hisat2_dir}/${SampleID}.hg38.sort.bam \
1>${Hisat2_dir}/${SampleID}.hisat2.log 2>&1


#get time end the job
echo "Job finished at:" `date`
