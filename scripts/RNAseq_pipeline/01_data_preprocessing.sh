#!/bin/bash

set -e
echo "Job Start at `date`"

##############Part1 QC: Rawdata to CleanData###########################

#Please copy the raw sequencing data into the directory "01_RawData"
mkdir 01_RawData 02_CleanData

#The path of the software
fastp=~/02.Software/miniconda3/bin/fastp
fastqc=~/02.Software/miniconda3/bin/fastqc

RawData_dir=./01_RawData
CleanData_dir=./02_CleanData
SampleID=SampleID

#filter low quality sequence by fastp 
${fastp} -i ${RawData_dir}/${SampleID}.fastq.gz \
-o ${CleanData_dir}/${SampleID}.clean.fastq.gz \
-h ${CleanData_dir}/${SampleID}.fastp.html \
-j ${CleanData_dir}/${SampleID}.fastp.json \
-w 4

#quality control by fastqc 
cd ${CleanData_dir}
${fastqc} -f fastq ${SampleID}.clean.fastq.gz

#return to initial directory
cd ..
#get time end the job
echo "Job finished at:" `date`