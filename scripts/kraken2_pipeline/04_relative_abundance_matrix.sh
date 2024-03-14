#!/bin/bash

set -e
echo "Job Start at `date`"

##############Part4 Relative abundance matrix###########################

mkdir 05_Relative_abundance

#The path of the software
combine_mpa=~/02.Software/KrakenTools-1.2/combine_mpa.py
relab=~/02.Software/miniconda3/bin/humann_renorm_table

kraken2_dir=./04_kraken2
relab_dir=./05_Relative_abundance

#combine_mpa
python ${combine_mpa} -i ${kraken2_dir}/*bracken.mpa \
-o ${relab_dir}/allsample.mpa

#species count table
cat ${relab_dir}/allsample.mpa | \
sed s/.bracken_report//g |  \
grep -E "s__|Classification" | \
grep -v Homo | \
sed 's/^.*s__//g' | \
sed s/#Classification/Classification/ > ${relab_dir}/species_readcount.txt

#genus count table
cat ${relab_dir}/allsample.mpa |  \
sed s/.bracken_report//g |  \
grep -E "g__|Classification" |  \
grep -v Homo |  \
grep -v "s__" |  \
sed 's/^.*g__//g' |  \
sed s/#Classification/Classification/ > ${relab_dir}/genus_readcount.txt

#species relative abundance
$relab -i ${relab_dir}/species_readcount.txt -u relab -o ${relab_dir}/species_abundances.txt

#genus species_abundances
$relab -i ${relab_dir}/genus_readcount.txt -u relab -o ${relab_dir}/genus_abundances.txt


#get time end the job
echo "Job finished at:" `date`