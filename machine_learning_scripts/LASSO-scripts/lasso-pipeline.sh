#!/bin/bash

ff=Host-Infection_vs_Cancer.diff.txt
train=trainset.2
test=testset.2
group=Infection

# repeat 100 
for i in {1..100}

do 
	Rscript mlr3_LASSO_cv_glmnet.R $ff $train $test $group $i
done


cat *selectfeature.txt | grep -v Intercept > allselect
Rscript feature_time.R
mkdir feature100
mv *_lasso_* feature100/

Rscript 2.mlr3_Lasso_model_SFS.R $ff $train $test $group
Rscript plotgenetime.R
