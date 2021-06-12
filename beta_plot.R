### Beta多样性作图，需要抽平后的OTU表
library("phyloseq")
library("ggsci")
library("vegan")
library("ggplot2")


setwd("D:/Users/lina1/Desktop/mfl")
#读取标准化后的OTU表,读取元数据(从qiime导出)
table <- read.table("Rare_12000_ASV.tsv",sep = "\t",header = T,check.names = FALSE)
rownames(table)<-table$Sample
table <-table[,-1]
metadata<-read.table("sample-metadata.tsv", header=T,row.names=1, sep="\t", comment.char="", stringsAsFactors=F)

main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=12),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=7),
                   text=element_text(family="sans", size=12))

##不同肠段排序图
table<-t(table)
physeq<-otu_table(table,taxa_are_rows = TRUE)
metadata_physeq<-sample_data(metadata)
phyloseq<-phyloseq(physeq,metadata_physeq)
ordinate <-ordinate(phyloseq,distance = "bray",method = "NMDS")

plot<-as.data.frame(cbind(ordinate[["points"]]))
tmp<- as.data.frame(cbind(phyloseq@sam_data[["sample"]],as.character(phyloseq@sam_data[["body.site"]])))
plot <- cbind(plot,tmp)
names(plot) <- c("NMDS1","NMDS2","Sample","region")
plot$region <- factor(plot$region,levels = c("ileum","cecum","colon","rectum"))

p<-ggplot(plot,aes(NMDS1,NMDS2,color=region,shape=region))+
  geom_point(size=2.6)+
  scale_shape_manual(values = c(17,16,16,16))+
  scale_color_manual(values = c("#4DBBD5E5","#E57272FF","#00A087E5","#7876B1B2","#feb24c"))+
  theme_bw()+
  main_theme

ggsave(file="GutREGION_NMDS.pdf",p,width=8,height=5)

stress <- ordinate$stress

##注意在AI中手动把stress加进去
