library("tidyr")
library("ggplot2")
library("ggsci")
library("reshape2")
setwd("D:/Users/lina1/Desktop/mfl/各水平OTU表")

#####属水平的堆积图#######
table <- read.delim("genus_table.txt")
 

#统计在每个样本中top10的属
top10_genus <- order(rowSums(table[,-1]),decreasing = T)[1:10]

#将在genus_top中的保留下来，其余的标记为others_table，并求和记为others,然后合并sub_table和others
sub_table=table[top10_genus,]
rownames(sub_table)<-sub_table[,1]
sub_table <- sub_table[,-1]
others_table=table[!rownames(table) %in% top10_genus,]
others=t(as.data.frame(colSums(others_table[,-1])))
table=rbind(sub_table,others)
rownames(table)[nrow(table)]<-"others"
rm(list = c("top10_genus","others","others_table","sub_table"))

#竖着排列
table<-cbind(rownames(table),table)
names(table)[1]<-"names"
plot=melt(table,id="names")
plot[,"region"] <- "ileum"
plot[,"region"][331:660] <- "colon"
plot[,"region"][661:990] <- "cecum"
plot[,"region"][991:1320] <- "rectum"
names(plot)=c('taxonomy',"sample","value","region")


#11个颜色配色
mypal<-c( "#E64B35FF","#0072B5FF","#4DBBD5FF","#20854EFF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
          "#91D1C2FF","#DC0000FF","#7E6148FF") 

#按照others大小对sample进行排序
levels <- subset(plot,taxonomy=="others",select = c(taxonomy,sample,value))
levels <- levels[order(levels$value,decreasing = F),]#others从小到大重排序
levels <- factor(levels$sample,levels = as.character(levels$sample))
plot$sample <- factor(plot$sample,levels = levels)

taxonomy_levels <- names(sort(rowSums(table[,-1]),decreasing = T))
taxonomy_levels <- c(taxonomy_levels[-1],"others")
plot$taxonomy <- factor(plot$taxonomy,levels = taxonomy_levels)
plot$region <- factor(plot$region,levels = c("ileum","cecum","colon","rectum"))
#可视化

theme <-theme(axis.text.x=element_blank(),
              axis.text.y=element_text(size=12),
              legend.text= element_text(size=10),
              title = element_text(size = 15),
              axis.ticks=element_blank()
              )


p1<-ggplot(plot,aes(x=sample, fill=taxonomy, y=value*100))+
  geom_col(position='stack',alpha=0.85,width = 0.9)+ #width=1 去掉间距
  scale_fill_manual(values = mypal)+
  labs(x='', y='Relative Abundance (%)')+
  scale_y_continuous(expand=c(0, 0))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+theme

p2<-p1+facet_wrap(~region,scales = "free_x",nrow = 2)

ggsave(p2,filename = "stackPlot.pdf",width = 12,height = 5)


