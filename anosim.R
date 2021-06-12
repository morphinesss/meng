library("vegan")

###利用adonis解释样品的差异度
setwd("D:/Users/lina1/Desktop/mfl")
#读取标准化后的OTU表,读取元数据(从qiime导出)
table <- read.table("Rare_12000_ASV.tsv",sep = "\t",header = T,check.names = FALSE)
rownames(table)<-table$Sample
table<-table[,-1]

##给数据分组
table<-cbind(table,rownames(table))
table <- separate(table,col = 'rownames(table)',into = c("tmp","region"))
table[table$region %in% seq(1,30),"region"] <- "ileum"
table[table$region %in% seq(31,60),"region"] <- "colon"
table[table$region %in% seq(61,90),"region"] <- "cecum"
table[table$region %in% seq(91,120),"region"] <- "rectum"
table$region<-factor(table$region,levels = c("ileum","cecum","colon","rectum"))

anosim<-anosim(subset(table,select = c(-tmp,-region)),table$region,permutations = 1000,distance = "bray")

##把R和p手动加到NMDS中





