#!/bin/bash

set -e
echo "Job Start at `date`"

##############Part4 expression profiles: Merge featurecount result of each sample into an expression profiles###########################

mkdir 05_Expression_profile

#The path of the script
merge=./merge_featurecount.py

featureCounts_dir=./04_featureCounts
Expression_profile_dir=./05_Expression_profile


#Expression_profile
python ${merge} -d ${featureCounts_dir} -o ${Expression_profile_dir}/RNAseq_featurecount.txt


#get time end the job
echo "Job finished at:" `date`
