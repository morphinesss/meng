###venn
library(tidyr)
library(pheatmap)
library(dplyr)

setwd("D:/Users/lina1/Desktop/mfl")

#导入数据
genus_table <- read.table("属水平.txt",header = T,sep = "\t")
lefse <- read.table("B_gutReion_lefse.txt",header = T,sep = "\t")
genus_table <- subset(genus_table,index %in% lefse$Biomark)
levels <- lefse$Biomark

heatmap <- NULL
list=c("ileum","colon","cecum","rectum")

number =c(30,60,90,120)
for (i in number){
  group<-subset(genus_table,select=seq(i-28,i+1))
  heatmap<-cbind(heatmap,assign(paste0(i,"_sum"),rowSums(group)/ncol(group)))
  rm(list=paste0(i,"_sum"))
}
heatmap <- as.data.frame(heatmap)
colnames(heatmap)<-list
heatmap <- heatmap[,c(1,3,2,4)]

heatmap <- as.data.frame(cbind(as.character(genus_table$index),heatmap))
names(heatmap)[1] <- "genus"
rownames(heatmap) <- heatmap$genus
heatmap <- heatmap[,-1]
heatmap[,c(1,2,3,4)] <- lapply(heatmap[,c(1,2,3,4)],as.character)
heatmap[,c(1,2,3,4)] <- lapply(heatmap[,c(1,2,3,4)],as.numeric)
heatmap <- heatmap[match(lefse[,1],rownames(heatmap)),]


pheatmap(heatmap,cluster_cols=F,cluster_rows=F,scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize_row = 6,
         main = "Heatmap showing gut regions associated genus",fontsize = 8,
         angle_col = 45,border_color = NA)


