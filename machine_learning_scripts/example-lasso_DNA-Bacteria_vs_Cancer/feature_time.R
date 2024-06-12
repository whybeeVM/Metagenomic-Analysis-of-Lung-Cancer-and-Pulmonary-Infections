library(ggpubr)
library(ggplot2)

data <- read.table("allselect", stringsAsFactors = F)
ti <- as.data.frame(table(data$V1))

ti <- ti[order(ti$Freq,decreasing = T),]
summary(ti$Freq)

ti$Var1 <- factor(ti$Var1, levels = rev(ti$Var1))

ggbarplot(ti, "Var1", "Freq",fill="#0073C2FF",
          color="#0073C2FF",
          palette = "Paired",label = T,
          rotate=T, xlab = "", lab.size = 3,
          lab.hjust = -0.1, lab.vjust = 0.5) +
  theme(
    axis.text = element_text(size = 9,colour="black"),
    legend.position = "none",
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
  ) + scale_y_continuous(limits = c(0,max(ti$Freq)+10), expand=c(0,0))

ggsave("feature.png", width = 7, height = 7)          
ggsave("feature.pdf", width = 7, height = 7)

write.table(ti, "freqorder.txt", sep = "\t", quote = F,
            row.names = F)
