library(future)
library(mlr3)
library(mlr3viz)
library(mlr3learners)
library(mlr3fselect)
library(mlr3filters)
library(ranger)
library(precrec)
library(ggplot2)
library(ggpubr)
library(pROC)
library(ROCR)
library(caTools)

args <- commandArgs(trailingOnly=TRUE)

file = args[1]
group_file = args[2]
group1 = args[4]
group2 = "Cancer"

testfile = args[1]
tegroup = args[3]

name = paste0(group1, "_vs_", group2, "_lasso")

#read data and processing
data <- read.table(file, row.names = 1, header = TRUE, sep  ="\t",
                   check.names = F, stringsAsFactors = F)
group <- read.table(group_file, header = F, check.names = F, 
                    stringsAsFactors = F)
colnames(group) <- c("sample", "groups")
rownames(group) <- group$sample

group <- subset(group, group$groups==group1 | group$groups==group2)

#diff <- read.table(diff_file, header = F, sep="\t",
#                   check.names = F, stringsAsFactors = F)
#rownames(diff) <- diff$V1
oversam <- intersect(rownames(group), colnames(data))
overfeature <- rownames(data)
sub_group <- group[oversam,]
sub_data <- as.data.frame(t(data[overfeature, oversam]))

train_set <- cbind(sub_data, sub_group)
train_set <- subset(train_set, select=-c(sample))

colnames(train_set) <- gsub("-","_", colnames(train_set))
colnames(train_set) <- gsub(" ","_", colnames(train_set))
colnames(train_set) <- gsub(":","_", colnames(train_set))

#test
testdata <- read.table(testfile, header = T,stringsAsFactors = F,sep="\t",
                       row.names = 1, check.names = F)
testdata <- testdata[overfeature,]
testgroup <- read.table(tegroup,row.names = 1,stringsAsFactors = F)
colnames(testgroup) <- "groups"
testgroup$sample <- rownames(testgroup)

testgroup <- subset(testgroup, testgroup$groups==group1 | testgroup$groups==group2)

testdata <- as.data.frame(t(testdata))
overlap <- intersect(rownames(testgroup), rownames(testdata))
sub_group <- testgroup[overlap,]
sub_otu <- testdata[overlap,]
test_set <- cbind(sub_otu, sub_group)
test_set <- subset(test_set, select=-c(sample))
colnames(test_set) <- gsub("-","_", colnames(test_set))
colnames(test_set) <- gsub(" ","_", colnames(test_set))
colnames(test_set) <- gsub(":","_", colnames(test_set))

#makeClassif train and feature select
set.seed(12345)
future::plan("multisession",workers = 10)

classif.lrn = lrn("classif.cv_glmnet", predict_type = "prob")

measure = list(msr("classif.auc"),msr("classif.acc"),
               msr("classif.precision"), msr("classif.recall"),
               msr("classif.sensitivity"), msr("classif.specificity"),
               msr("classif.tpr"), msr("classif.fpr"))

featureorder <- read.table("freqorder.txt", header = T, stringsAsFactors = F)

allauc <- data.frame()
for (i in 2:nrow(featureorder)) {
  select_ft <- featureorder[1:i,1]
  fs_table <- train_set[,c(select_ft, "groups")]
  fs.task = as_task_classif(fs_table, target = "groups", id="groups")
  fs.task$col_roles$stratum = fs.task$col_roles$target
  fs.task
  
  #feature select
  classif.lrn$train(fs.task)
  classif.lrn$model
  
  pdf(paste0(name, "_", i, "_LASSO_cv_lambda.pdf"), height = 5, width = 6)
  plot(classif.lrn$model)
  dev.off()
  
  png(paste0(name, "_", i,"_LASSO_cv_lambda.png"), height = 5, width = 6,units="in", res=350)
  plot(classif.lrn$model)
  dev.off()
  
  lambda_parm = classif.lrn$model$lambda.min
  feature_coef = coef(classif.lrn$model, s = "lambda.min") 
  feature_coef
  select_feature = feature_coef@Dimnames[[1]][feature_coef@i+1]
  select_feature
  
  write.table(select_feature, paste0(name, "_", i,"_selectfeature.txt"),
              sep = "\t", row.names = F, quote = F, col.names = F)
  
  #final model
  classif.lrn$param_set$values$s = lambda_parm
  pred = classif.lrn$predict(fs.task)
  pred
  pred$confusion
  
  pred_res = data.frame(row_ids=pred$row_ids,
                        truth=pred$truth, response=pred$response,
                        prob.positive=pred$prob[,1],
                        prob.negative=pred$prob[,2])
  
  write.table(pred_res, paste0(name, "_", i,"_Trainset_pred.txt"),
              sep = "\t", row.names = F, quote = F)
  
  # predict test set
  test_fs = test_set[,c(select_ft, "groups")]
  test.task = as_task_classif(test_fs, target = "groups", id="groups")
  test.task
  test_result = classif.lrn$train(fs.task)$predict(test.task)
  test_pred_res = data.frame(row_ids=test_result$row_ids,
                             truth=test_result$truth, 
                             response=test_result$response,
                             prob.positive=test_result$prob[,1],
                             prob.negative=test_result$prob[,2])
  
  write.table(test_pred_res, paste0(name, "_", i,"_Testset_pred.txt"),
              sep = "\t", row.names = F, quote = F)
  
  #plot AUC 
  rocobj <- roc(pred_res$truth,pred_res[,4],auc = TRUE,
                ci=TRUE, print.auc=TRUE,direction=">",levels=c(fs.task$positive, fs.task$negative))
  rocobj_test <- roc(test_pred_res$truth,test_pred_res[,4],auc = TRUE,
                     ci=TRUE, print.auc=TRUE,direction=">",levels=c(test.task$positive, test.task$negative))
  ciobj <- ci.se(rocobj, specificities=seq(0, 1, 0.01))
  auc <- auc(rocobj)[1]
  auc_low <- ci(rocobj,of="auc")[1]
  auc_high <- ci(rocobj,of="auc")[3]
  auc_full <- paste0("Train set AUC=",round(auc,digits = 3)," (",
                     round(auc_low,digits = 3),"-",round(auc_high,digits = 3),")")
  auc_full
  auc_print <- paste0(group1, "-", group2, "trainset AUC = ",round(auc,digits = 3),"(95% CI:",
                      round(auc_low,digits = 3),"-",round(auc_high,digits = 3),")", " testset AUC=", rocobj_test$auc)
  
  write.table(auc_print, paste0(name, "_", i, "_AUC.txt"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  data_ci <- ciobj[1:101,1:3]
  data_ci <- as.data.frame(data_ci)
  x = as.numeric(rownames(data_ci))
  data_ci <- data.frame(x,data_ci)
  
  trainplot <- data.frame(rocobj$specificities, rocobj$sensitivities)
  colnames(trainplot) <- c("Specificity","Sensitivity")
  trainplot$method <- auc_full
  testplot <- data.frame(rocobj_test$specificities, rocobj_test$sensitivities)
  colnames(testplot) <- c("Specificity","Sensitivity")
  testplot$method <- paste0("Validation set AUC = ",round(rocobj_test$auc, 3))
  
  classif_result <- rbind(trainplot, testplot)
  #plot AUC
  
  p <- ggplot(classif_result) + 
    geom_path(aes(x = 1-Specificity, y = Sensitivity, color = method), linewidth=1) +
    theme_bw() +
    #scale_x_reverse()+
    scale_color_manual(values = c("#E64B35FF","#4DBBD5FF")) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 color="darkgrey", linetype=6)+
    geom_ribbon(data = data_ci,aes(x=1-x,ymin=X2.5.,ymax=X97.5.),
                fill = 'grey',alpha=0.3)+
    labs(x="1-Specificity", y="Sensitivity", title = "ROC curve", color = "") +
    theme(plot.title=element_text(hjust=0.5),
          legend.position = c(0.68,0.15),
          legend.background = element_rect(fill=NA),
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
          axis.text.x = element_text(color = "black",size = 12),
          axis.text.y = element_text(color = "black",size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14))
  
  ggsave(paste0(name, "_", i,"_Train_Test_AUC.png"), p, width = 5, height = 5)
  ggsave(paste0(name, "_", i,"_Train_Test_AUC.pdf"), p, width = 5, height = 5)
  
  tmp <- as.data.frame(c(i,auc(rocobj)[1],ci(rocobj,of="auc")[1],ci(rocobj,of="auc")[3],
           auc(rocobj_test)[1],ci(rocobj_test,of="auc")[1],ci(rocobj_test,of="auc")[3],
           paste0(select_ft, collapse = "|")))
  rownames(tmp) <- c("number","trainAUC","trainCI1","trainCI2",
                     "testAUC","testCI1","testCI2","feature")
  colnames(tmp) <- i
  allauc <- rbind(allauc, t(tmp))
}
write.table(allauc, "allfeatureAUC.txt", sep = "\t", row.names = F,
            quote = F)
write.csv(allauc, "allfeatureAUC.csv", row.names = F, quote = F)
