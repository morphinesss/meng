library("vegan")
library("ggplot2")
library("ggsci")
library("agricolae")
library("ggpubr")
setwd("D:/Users/lina1/Desktop/mfl")

#读取标准化后的OTU表
table <- read.table("Rare_12000_ASV.tsv",sep = "\t",header = T,check.names = FALSE)
rownames(table)<-table$Sample
table<-table[,-1]

#定义计算多样性指数的函数
alpha_index <- function(x,method = 'richness',tree = NULL,base=exp(1)){ #对数底数使用e
  if (method == "richness") result <- rowSums(x>0)  #丰富度指数
  else if (method == "chao1") result <- estimateR(x)[2, ]    #Chao1 指数
  else if (method == 'ace') result <- estimateR(x)[4, ]    #ACE 指数
  else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)
  else if (method == 'simpson') result <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou 均匀度
  result  #函数将会返回最后一行运行的结果
}

#计算多样性指数
shannon <- alpha_index(table,method = "shannon",tree=NULL,base=exp(1))
richness <- alpha_index(table,method = "richness",tree=NULL,base=exp(1))
chao1 <- alpha_index(table,method = "chao1",tree=NULL)
simpson <- alpha_index(table,method = "simpson")

sample<-names(shannon)
alpha_index_table <- as.data.frame(cbind(sample,shannon,richness,chao1,simpson))

#写出α多样性指数
#write.table(file="alpha_index_table.txt",alpha_index_table,sep="\t",quote = FALSE,row.names = FALSE)

alpha_index_table<-cbind(alpha_index_table,alpha_index_table$sample)
alpha_index_table <- separate(alpha_index_table,col = 'alpha_index_table$sample',into = c("tmp","region"))
alpha_index_table[alpha_index_table$region %in% seq(1,30),"region"] <- "ileum"
alpha_index_table[alpha_index_table$region %in% seq(31,60),"region"] <- "colon"
alpha_index_table[alpha_index_table$region %in% seq(61,90),"region"] <- "cecum"
alpha_index_table[alpha_index_table$region %in% seq(91,120),"region"] <- "rectum"
alpha_index_table$region<-factor(alpha_index_table$region,levels = c("ileum","cecum","colon","rectum"))
alpha_index_table <- alpha_index_table[,-6]



#先将因子型转换为字符型在转换成数值型
alpha_index_table$shannon<-as.numeric(as.character(alpha_index_table$shannon)) 
alpha_index_table$richness<-as.numeric(as.character(alpha_index_table$richness))
alpha_index_table$chao1<-as.numeric(as.character(alpha_index_table$chao1))
alpha_index_table$simpson<-as.numeric(as.character(alpha_index_table$simpson))

#画图
mypal<-c("#00468BE5","#B71B1BFF","#20854EB2","#A20056FF")

theme <- theme(legend.text= element_text(size=13),
               axis.ticks=element_line(color="black"),
               axis.text=element_text(color="black", size=13),
               axis.line.x=element_line(size=.5, colour="black"),
               axis.line.y=element_line(size=.5, colour="black"),
               axis.title = element_text(face = "bold",size = 14),
               axis.text.x = element_text(angle=45, hjust=1,size = 13),
               axis.text.y = element_text(size = 13),
               legend.position="right")

p2<-ggboxplot(alpha_index_table,x="region",y="shannon",color="region",
              palette = mypal,add = "jitter",xlab = "",ylab = "Shannon index")+
  theme_bw()+theme
p
##stat_compare_means() 添加星号显著性标记
mycompare=list(c("ileum","cecum"),c("ileum","colon"),c("ileum","rectum"),c("cecum","colon"),c("cecum","rectum"),c("colon","rectum"))
p2<-p2+stat_compare_means(comparisons = mycompare,label = "p.signif",method = "wilcox")
ggsave(p2,file="shannon.pdf",width=5,height=5)



  

