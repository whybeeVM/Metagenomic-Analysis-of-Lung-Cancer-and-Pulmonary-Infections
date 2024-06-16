library(ggplot2)
library(ggpubr)
library(pROC)
library(ROCR)

args <- commandArgs(trailingOnly=TRUE)

allfile <- read.table("file.txt", stringsAsFactors = F)

allauc <- data.frame()
# CNV DNA RNA Host TE
tt <- "all"

model <- c("Lasso","RF","SVM","XGBoost")
for (mm in 1:4) {
  for (nn in 2:4) {
    if (nn > mm) {
      n1 <- model[mm]
      n2 <- model[nn]
      
      #dir <- paste0(tt, "_combined_",n1,"_",n2)
      #dir.create(dir)
      
      ff <- allfile[c(grep(n1, allfile$V1), grep(n2, allfile$V1)),]
      ff <- ff[grepl(tt,ff)]
      
      for (name in c("Bacteria_vs_Cancer","Fungi_vs_Cancer",
                     "TB_vs_Cancer","Infection_vs_Cancer")) {
        
        subff <- ff[grepl(name, ff)]
        
        model1 <- read.table(subff[2], header = T, stringsAsFactors = F)
        model2 <- read.table(subff[4], header = T, stringsAsFactors = F)
        
        model1_test <- read.table(subff[1], header = T, stringsAsFactors = F)
        model2_test <- read.table(subff[3], header = T, stringsAsFactors = F)
        
        #train prob
        model1_mean <- aggregate(model1$prob.positive, list(row_ids=model1$row_ids,truth=model1$truth), mean)
        colnames(model1_mean)[3] <- "prob.positive"
        model1_mean$prob.negative <- 1 - model1_mean$prob.positive
        
        model2_mean <- aggregate(model2$prob.positive, list(row_ids=model2$row_ids,truth=model2$truth), mean)
        colnames(model2_mean)[3] <- "prob.positive"
        model2_mean$prob.negative <- 1 - model2_mean$prob.positive
        
        model1_mean <- model1_mean[order(model1_mean$row_ids), ]
        model2_mean <- model2_mean[order(model2_mean$row_ids), ]
        
        #test prob
        model1_test <- model1_test[order(model1_test$row_ids), ]
        model2_test <- model2_test[order(model2_test$row_ids), ]
        
        #combind
        
        for (i in seq(0,1,0.1)) {
          
          combind <- data.frame("row_ids"=model2_mean$row_ids, "truth"=model2_mean$truth)
          combind$posi <- i * model1_mean$prob.positive + (1-i) * model2_mean$prob.positive
          combind$nega <- i * model1_mean$prob.negative + (1-i) * model2_mean$prob.negative
          
          combind_test <- data.frame("row_ids"=model2_test$row_ids, "truth"=model2_test$truth)
          combind_test$posi <- i * model1_test$prob.positive + (1-i) * model2_test$prob.positive
          combind_test$nega <- i * model1_test$prob.negative + (1-i) * model2_test$prob.negative
          
          # write.table(combind, paste0(dir ,"/", name, "_",i,"_Train_pred.txt"), sep = "\t", row.names = F,
          #             quote = F)
          # write.table(combind_test, paste0(dir ,"/", name, "_",i,"_Test_pred.txt"), sep = "\t", row.names = F,
          #             quote = F)
          # 
          g2 <- gsub("_vs_.*","", name)
          g1 <- gsub(".*_vs_", "", name)
          if(g1 == sort(c(g1,g2))[1]) {
            direct = "<"
          }else{
            direct = ">"
          }
          
          rocobj <- roc(combind$truth,combind$posi,auc = TRUE,
                        ci=TRUE, print.auc=TRUE,direction = direct,
                        levels=c(g2,g1))
          rocobj_test <- roc(combind_test$truth,combind_test$posi,auc = TRUE,
                             ci=TRUE, print.auc=TRUE,direction = direct,
                             levels=c(g2,g1))
          if (rocobj$auc < 0.5) {
            rocobj <- roc(combind$truth,combind$nega,auc = TRUE,
                          ci=TRUE, print.auc=TRUE, direction = direct,
                          levels=c(g2,g1))
          }
          if (rocobj_test$auc<0.5) {
            rocobj_test <- roc(combind_test$truth,combind_test$nega,auc = TRUE,
                               ci=TRUE, print.auc=TRUE, direction = direct,
                               levels=c(g2,g1))
          }
          
          roc_result <- coords(rocobj, "best",
                               ret=c("threshold", "specificity", 
                                     "sensitivity","accuracy"))
          roctest_result <- coords(rocobj_test, "best",
                               ret=c("threshold", "specificity", 
                                     "sensitivity","accuracy"))
          if (nrow(roc_result) > 1) {
            roc_result <- coords(rocobj, 0.5,
                                 ret=c("threshold", "specificity", 
                                       "sensitivity","accuracy"))
          }
          if(nrow(roctest_result) > 1) {
            roctest_result <- coords(rocobj_test, 0.5,
                                     ret=c("threshold", "specificity", 
                                           "sensitivity","accuracy"))
            
          }
          
        #   ciobj <- ci.se(rocobj, specificities=seq(0, 1, 0.01))
        #   auc <- auc(rocobj)[1]
        #   auc_low <- ci(rocobj,of="auc")[1]
        #   auc_high <- ci(rocobj,of="auc")[3]
        #   auc_full <- paste0("Train set AUC=",round(auc,digits = 3)," (",
        #                      round(auc_low,digits = 3),"-",round(auc_high,digits = 3),")")
        #   auc_full
        # 
        #   data_ci <- ciobj[1:101,1:3]
        #   data_ci <- as.data.frame(data_ci)
        #   x = as.numeric(rownames(data_ci))
        #   data_ci <- data.frame(x,data_ci)
        #   
        #   trainplot <- data.frame(rocobj$specificities, rocobj$sensitivities)
        #   colnames(trainplot) <- c("Specificity","Sensitivity")
        #   trainplot$method <- auc_full
        #   testplot <- data.frame(rocobj_test$specificities, rocobj_test$sensitivities)
        #   colnames(testplot) <- c("Specificity","Sensitivity")
        #   testplot$method <- paste0("Validation set AUC = ",round(rocobj_test$auc, 3))
        #   
        #   classif_result <- rbind(trainplot, testplot)
        #   #plot AUC
        #   
        #   p <- ggplot(classif_result) + 
        #     geom_path(aes(x = 1-Specificity, y = Sensitivity, color = method), linewidth=1) +
        #     theme_bw() +
        # #    scale_x_reverse()+
        #     scale_color_manual(values = c("#E64B35FF","#4DBBD5FF")) +
        #     geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
        #                  color="darkgrey", linetype=6)+
        #     geom_ribbon(data = data_ci,aes(x=1-x,ymin=X2.5.,ymax=X97.5.),
        #                 fill = 'grey',alpha=0.3)+
        #     labs(x="1-Specificity", y="Sensitivity", title = "ROC curve", color = "") +
        #     theme(plot.title=element_text(hjust=0.5),
        #           legend.position = c(0.68,0.15),
        #           legend.background = element_rect(fill=NA),
        #           panel.grid.major =element_blank(),
        #           panel.grid.minor = element_blank(),
        #           plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        #           axis.text.x = element_text(color = "black",size = 12),
        #           axis.text.y = element_text(color = "black",size = 12),
        #           axis.title.x = element_text(size = 14),
        #           axis.title.y = element_text(size = 14))
        #   p
        #   ggsave(paste0(dir, "/", name, "_",i, "_Train_Test_AUC.png"), p, width = 5, height = 5)
        #   ggsave(paste0(dir, "/", name, "_",i, "_Train_Test_AUC.pdf"), p, width = 5, height = 5)
          if (i == 0) {
            modelname <- n2 
          }else if( i == 1){
            modelname <- n1
          }else{
            modelname <- paste0(i, "*", n1, "+", (1-i), "*", n2)
          }
          tmp <- as.data.frame(c(modelname, name,
                   auc(rocobj)[1],ci(rocobj,of="auc")[1],ci(rocobj,of="auc")[3],
                   roc_result$threshold, roc_result$sensitivity, roc_result$specificity,
                   roc_result$accuracy,
                   auc(rocobj_test)[1],ci(rocobj_test,of="auc")[1],ci(rocobj_test,of="auc")[3],
                   roctest_result$threshold, roctest_result$sensitivity,
                   roctest_result$specificity, roctest_result$accuracy))
          rownames(tmp) <- c("model","compare","trainAUC","trainCI1","trainCI2",
                             "train_threshold","train_sensitivity","train_specificity",
                             "train_accuracy",
                             "testAUC","testCI1","testCI2",
                             "test_threshold","test_sensitivity","test_specificity",
                             "test_accuracy")
          colnames(tmp) <- NULL
          allauc <- rbind(allauc, t(tmp))
        }
     
      }
    }
  }
}

write.table(allauc, paste0(tt, "_combinedAUC_SEN_SPE.txt"), sep = "\t", row.names = F,
            quote = F)
write.csv(allauc, paste0(tt, "_combinedAUC_SEN_SPE.csv"), row.names = F, quote = F)

#heatmap
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)
library(aplot)
library(ggtree)

allauc$test_sensitivity <- as.numeric(allauc$test_sensitivity)
allauc$test_specificity <- as.numeric(allauc$test_specificity)
allauc$test_accuracy <- as.numeric(allauc$test_accuracy)
allauc$testAUC <- as.numeric(allauc$testAUC)

# 14 15 16
colnum <- c(10,14,15,16)
colname <- c("AUC","Sensitivity","Specificity","Accuracy")

for (kkk in 1:4) {

  testauc <- allauc[,c(1,2,colnum[kkk])]
  titlename <- colname[kkk]
  
  colnames(testauc)[3] <- "testAUC"
  testauc <- testauc[!duplicated(testauc), ] 
  testauc_matix <- as.data.frame(acast(testauc, testauc$model~testauc$compare))
  testauc_mean <- as.data.frame(apply(testauc_matix,1, mean))
  colnames(testauc_mean) <- "MeanAUC"
  testauc_mean$name <- rownames(testauc_mean)
  testauc_mean <- as.data.frame(testauc_mean[order(testauc_mean$MeanAUC, decreasing = T),])
  
  #plot mean auc
  testauc_mean$name <- factor(testauc_mean$name, levels = rev(testauc_mean$name))
  testauc_mean$MeanAUC <- round(testauc_mean$MeanAUC,3)
  testauc_mean$col <- "black"
  testauc_mean[1,"col"] <- "red"
  
  p1 <- ggbarplot(testauc_mean, "name", "MeanAUC",fill="#0073C2FF",
            color="#0073C2FF",
            palette = "Paired",label = T, 
            title = paste0("Mean ",titlename),# lab.col = testauc_mean$col,
            rotate=T, xlab = "",ylab="", lab.size = 2.5,
            lab.hjust = -0.1, lab.vjust = 0.5) +
    theme_gray() +
    theme(
      axis.text = element_text(size = 9,colour="black"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
    ) + scale_y_continuous(limits = c(0,max(testauc_mean$MeanAUC)+0.3),
                           expand=c(0,0))
  p1
  
  #heatmap figure
  groupcolor <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF")
  order <- c("Infection_vs_Cancer","Bacteria_vs_Cancer",
             "Fungi_vs_Cancer","TB_vs_Cancer")
  
  testauc$model <- factor(testauc$model, levels = rev(rownames(testauc_mean)))
  testauc$compare <- factor(testauc$compare, levels = order)
  
  pheat <- ggplot(testauc, aes(compare,model,fill=testAUC))+ geom_tile()+
    geom_text(aes(x = compare, y = model, label = round(testAUC,3)),
              size=3) +
    scale_y_discrete(position = "right")+
    #scale_fill_gradientn(colors = colorRampPalette(colors = c('#21b6af','white','#eeba4d'))(100)) +
    scale_fill_gradientn(colors = colorRampPalette(colors = rev(brewer.pal(11,"RdBu")[3:9]))(100)) + 
    labs(x=NULL,y=NULL,fill=titlename)+
    theme_minimal() +
    theme(
      #axis.text.y = element_text(size = 9,colour=c(rep("black",57),"red")),
      axis.text.y = element_text(size = 9,colour="black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank())
  pheat
  
  pheat1 <- pheat + scale_fill_gradientn(colors = colorRampPalette(colors = c('#21b6af','white','#eeba4d'))(100)) 
  pheat1
  
  df_anno <-  data.frame(compare = order,
                         col = c(1,1,1,1))
  df_anno$compare <- factor(df_anno$compare, levels = order)
  pheat_anno <- ggplot(df_anno,aes(compare,col,fill=compare))+geom_tile()+
    scale_y_discrete(position = "right")+
    labs(x=NULL,y=NULL)+
    theme_minimal()+
    scale_fill_manual(values = groupcolor) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
          )
  pheat_anno
  
  hc <- hclust(dist(testauc_matix,method="maximum")) # 
  
  # phr <- ggtree(hc) + geom_text(aes(label=node)) + geom_tiplab(align = T) 
  # x_range <- 4 * layer_scales(phr)$x$range$range
  # phr <- phr + xlim(c(NA, x_range[2]))
  # ggsave(paste0(tt,"_tree.pdf"), phr, width = 10 ,height = 10)
  
  phr2 <- ggtree(hc) #+ geom_text(aes(label=node), size=2.5)
  # phr2 <- ggtree::rotate(phr2,69) # all kkk 3
  # phr2 <- ggtree::rotate(phr2,89)
  # phr2 <- ggtree::rotate(phr2,73)
  # phr2 <- ggtree::rotate(phr2,76)
  # phr2 <- ggtree::rotate(phr2,62)
  # phr2 <- ggtree::rotate(phr2,71)
  
  # phr2 <- ggtree::rotate(phr2,69)  # all kkk 1
  # phr2 <- ggtree::rotate(phr2,76)
  # phr2 <- ggtree::rotate(phr2,64)
  # phr2 <- ggtree::rotate(phr2,82)
  
  # phr2 <- ggtree::rotate(phr2,59)
  # phr2 <- ggtree::rotate(phr2,61)
  # phr2 <- ggtree::rotate(phr2,75)
  # phr2 <- ggtree::rotate(phr2,60)
  # phr2 <- ggtree::rotate(phr2,68)
  # 
  pp <- pheat %>%  
    insert_left(phr2,width = 0.2) %>%
    insert_top(pheat_anno,height = 0.04) %>% 
    insert_right(p1,width = 0.5)
  pp
  ggsave(paste0(tt, "_combined",titlename,"_heatmap1.png"), pp, width = 10, height = 10)
  ggsave(paste0(tt, "_combined",titlename,"_heatmap1.pdf"), pp, width = 10, height = 10)
  
  pp1 <- pheat1 %>%  
    insert_left(phr2,width = 0.2) %>%
    insert_top(pheat_anno,height = 0.04) %>% 
    insert_right(p1,width = 0.5)
  
  ggsave(paste0(tt, "_combined",titlename,"_heatmap2.png"), pp1, width = 10, height = 10)
  ggsave(paste0(tt, "_combined",titlename,"_heatmap2.pdf"), pp1, width = 10, height = 10)

}
