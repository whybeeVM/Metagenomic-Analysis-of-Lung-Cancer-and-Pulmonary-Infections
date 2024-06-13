#!/bin/bash

set -e
echo "Job Start at `date`"

##############Part3 kraken2 and bracken###########################

mkdir 04_kraken2

#The path of the software
kraken2=~/02.Software/miniconda3/bin/kraken2
bracken=~/02.Software/miniconda3/bin/bracken
k2mpa=~/02.Software/KrakenTools-1.2/kreport2mpa.py

#kraken2 db
db=~/03.Database/k2_pluspf_20210517/

rmhost_dir=./03_rmhost
kraken2_dir=./04_kraken2
SampleID=SampleID

#kraken2
${kraken2} -db ${db} \
--threads 10 \
--report ${kraken2_dir}/${SampleID}.report \
--use-names ${rmhost_dir}/${SampleID}.rmhost.fastq.gz  \ 
> ${kraken2_dir}/${SampleID}.kraken2

#bracken
${bracken} -d ${db} \
-i ${kraken2_dir}/${SampleID}.report \
 -o ${kraken2_dir}/${SampleID}.bracken.Abundance.Estimation \
 -w ${kraken2_dir}/${SampleID}.bracken_report -r 50

 #kreport2mpa
 python ${k2mpa} -r ${kraken2_dir}/${SampleID}.bracken_report \
 -o ${kraken2_dir}/${SampleID}.bracken.mpa --display-header

#get time end the job
echo "Job finished at:" `date`
