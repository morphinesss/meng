#多样性分析，利用vegan计算多样性指数
library("vegan")
library("ggplot2")
library("doBy")
#(1)设置工作路径，读取数据,并转置
setwd("D:/Users/lina1/Desktop/mfl")
table <- as.data.frame(read.table("feature-table.txt",sep = "\t",header = T,check.names = F),stringsAsFactors = F)
rownames(table)<-table$OTU.ID
table <- table[,-1]
vegan_table=t(table)

#定义函数，计算多种Alpha多样性指数，结果返回至向量
#function创建函数，基本语法为function(#arguments){#do something}
alpha_index <- function(x,method = 'richness',tree = NULL,base=exp(1)){ #对数底数使用e
  if (method == "richness") result <- rowSums(x>0)  #丰富度指数
  else if (method == "chao1") result <- estimateR(x)[2, ]    #Chao1 指数
  else if (method == 'ace') result <- estimateR(x)[4, ]    #ACE 指数
  else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)
  else if (method == 'simpson') result <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou 均匀度
  result  #函数将会返回最后一行运行的结果
}

#根据抽样步长(step),统计每个稀释梯度下的Alpha多样性指数，结果返回至列表
alpha_curves <- function(x, step, method = 'richness', rare = NULL, tree = NULL, base = exp(1)) {
  x_nrow <- nrow(x) #行的数目
  if (is.null(rare)) rare <- rowSums(x) else rare <- rep(rare, x_nrow)  #计算每个样本的OTU总条数
  alpha_rare <- list()
  
  for (i in 1:x_nrow) {
    step_num <- seq(0, rare[i], step) #每个样本按照步数抽样  
    if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i]) #如果step_num中最大的数值小于样本的OTU数，在step_num后面加上OTU数目
    
    #计算每个梯度下的多样性指数
    alpha_rare_i <- NULL
    for (step_num_n in step_num) alpha_rare_i <- c(alpha_rare_i, alpha_index(x = rrarefy(x[i, ], step_num_n), method = method, tree = tree, base = base))
    #重抽样函数rrarefy,rrarefy(x[i, ], step_num_n)返回矩阵，是重抽样结果
    names(alpha_rare_i) <- step_num
    alpha_rare <- c(alpha_rare, list(alpha_rare_i))
  }
  
  names(alpha_rare) <- rownames(x)
  alpha_rare #返回alpha_rare列表，
}

##基于一次抽样的结果可能存在着误差，希望在同等深度下多次抽样并统计均值和标准差，最终绘制成带有误差棒的曲线图
plot_richness <- data.frame()
for (n in 1:5){
  richness_curves <- alpha_curves(vegan_table,step = 5000,method="richness")
  for(i in names(richness_curves)){
    richness_curves_i <- (richness_curves[[i]])
    richness_curves_i <- data.frame(rare=names(richness_curves_i), alpha = richness_curves_i, sample = i, stringsAsFactors = FALSE)
    plot_richness <- rbind(plot_richness, richness_curves_i)
  }
}  

##计算均值和标准差（summaryBy()函数）,此分组中样本太多，不使用误差棒
plot_richness_stat <- summaryBy(alpha~sample+rare, plot_richness, FUN = c(mean, sd))
plot_richness_stat$rare <- as.numeric(plot_richness_stat$rare)
plot_richness_stat[which(plot_richness_stat$rare == 0),'alpha.sd'] <- NA #标准差为0写为NA

##设定分组
plot_richness_stat<-cbind(plot_richness_stat,plot_richness_stat$sample)
plot_richness_stat <- separate(plot_richness_stat,col = 'plot_richness_stat$sample',into = c("tmp","region"))
plot_richness_stat[plot_richness_stat$region %in% seq(1,30),"region"] <- "ileum"
plot_richness_stat[plot_richness_stat$region %in% seq(31,60),"region"] <- "colon"
plot_richness_stat[plot_richness_stat$region %in% seq(61,90),"region"] <- "cecum"
plot_richness_stat[plot_richness_stat$region %in% seq(91,120),"region"] <- "rectum"
plot_richness_stat$region<-factor(plot_richness_stat$region,levels = c("ileum","cecum","colon","rectum"))
plot_richness_stat <- plot_richness_stat[,-5]

##利用ggplot2绘制稀疏曲线
mypal<-c("#00468BE5","#B71B1BFF","#20854EB2","#A20056FF")

p1<-ggplot(plot_richness_stat, aes(rare, alpha.mean,group=sample,color=region))+
  geom_smooth(size=0.3,se=FALSE,span=1.2,method = 'loess',formula='y ~ x')+
  labs(x = 'Number of sequences', y = 'Richness', color = NULL)+
  scale_color_manual(values=mypal)+
  theme(panel.grid = element_blank(),panel.background = element_rect(fill = 'transparent', color = 'black'),legend.key = element_rect(fill = 'transparent'))
p2<-p1+facet_wrap(.~region,nrow = 2) #facet_wrap进行分面
ggsave(p2,file="Rarefaction_curves.pdf",width=10,height=10)
write.table(plot_richness_stat,"Rarefaction_curves.txt",quote = FALSE,row.names = F,sep = "\t")

##根据稀疏曲线的结果，选择在25000这个序列深度进行标准化,得到标准化之后的OTU表格
##删除OTU总条数小于25000的样本
new_vegan_table <- vegan_table[-which(rowSums(vegan_table)<=12000),]
x_nrow <- nrow(new_vegan_table) #现有样本的数目
rare <- rowSums(new_vegan_table) #计算现有每个样本的OTU总条数
Sample <- names(rare)
rarefaction_25000 <- NULL
for(i in 1:x_nrow){
  rarefaction_i = rrarefy(new_vegan_table[i, ], 12000) #重抽样
  rarefaction_25000<-rbind(rarefaction_25000,rarefaction_i)
}
rarefaction_25000 <- cbind(Sample,rarefaction_25000)
#写出标准化后的OTU表
write.table(file="Rare_12000_ASV.tsv",rarefaction_25000,col.names = T,row.names = F,sep = "\t",quote = FALSE)












