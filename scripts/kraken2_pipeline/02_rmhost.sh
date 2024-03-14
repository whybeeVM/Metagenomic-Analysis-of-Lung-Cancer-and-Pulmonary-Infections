#!/bin/bash

set -e
echo "Job Start at `date`"

##############Part2 filt host sequence###########################

mkdir 03_rmhost

#The path of the software
bwa=~/02.Software/miniconda3/bin/bwa
samtools=~/02.Software/miniconda3/bin/samtools
seqkit=~/02.Software/miniconda3/bin/seqkit

hg38fasta=hg38.fa
CleanData_dir=./02_CleanData
rmhost_dir=./03_rmhost
SampleID=SampleID

##index of reference genome
${bwa} index ${hg38fasta}

#filt host sequence
${bwa} mem -t 20  \
hg38.fa ${CleanData_dir}/${SampleID}.clean.fastq.gz | \
${samtools} sort --threads 8 -m 4G | \
${samtools} view --threads 8 -b -h -f 4 \
-o ${rmhost_dir}/${SampleID}.bam \
1>${rmhost_dir}/${SampleID}.log 2>&1

#exact rmhost sequence
${samtools} fastq --threads 8 \
${rmhost_dir}/${SampleID}.bam > \
${rmhost_dir}/${SampleID}.rmhost.fastq

gzip ${rmhost_dir}/${SampleID}.rmhost.fastq

#stat
${seqkit} stat ${rmhost_dir}/${SampleID}.rmhost.fastq.gz > ${rmhost_dir}/${SampleID}.stat


#get time end the job
echo "Job finished at:" `date`