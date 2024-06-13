library(ggplot2)
library(ggpubr)
library(ggsci)

data <- read.table("allfeatureAUC.txt", header = T, check.names = F,
                   stringsAsFactors = F)
maxauc <- data[which.max(data$trainAUC),1]
maxauc

ggplot(data, aes(x= number, y = trainAUC)) +
  geom_rect(aes(xmin=maxauc-0.5,xmax=maxauc+0.5,ymin=-Inf,ymax=Inf),
            fill="grey80",alpha=0.3)+
  geom_errorbar(aes(ymin = trainCI1, ymax =trainCI2), 
                width = 0.2, position = position_dodge(0.1) ) +
  geom_line(position = position_dodge(0.1),color="red") +
  geom_point(position = position_dodge(0.1),color="red") +
  theme_bw() +
  theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black",size = 12,),
        axis.title = element_text(color = "black",size = 12),
        legend.position = 'none') + 
  labs(x="Feature Number", y="AUC of train set") 

ggsave("SFSmaxAuc.png", width = 7, height = 5)
ggsave("SFSmaxAuc.pdf", width = 7, height = 5)

system("mkdir select")
system(paste0("cp *", maxauc, "*", " select"))

