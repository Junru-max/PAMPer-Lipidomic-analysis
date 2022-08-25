### correlation matrix and network 
library(ComplexHeatmap)
library(scales)
library(circlize)
library(dendsort)
library(Seurat)
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)

##Fig1B,Fig4A
######
lipid_class<-readRDS(file = "rdata/lipid_class.rdata")
lipid_class_stat <- lipid_class %>% 
  group_by(Class) %>% 
  count() %>% 
  ungroup()  %>% 
  mutate(per=`n`/sum(`n`)) %>% 
  arrange(desc(n))
lipid_class_stat$Class<-as.character(lipid_class_stat$Class)
lipid_class_stat$Class<-factor(lipid_class_stat$Class,levels =lipid_class_stat$Class )
lipid_class_stat$Species<-paste(lipid_class_stat$Class,lipid_class_stat$n,sep = ":")
lipid_class_stat$Species<-factor(lipid_class_stat$Species,levels =lipid_class_stat$Species )
ggplot(data=lipid_class_stat)+
  geom_bar(aes(x="", y=n, fill=Species), stat="identity", width = 1)+
  coord_polar("y", start=0, direction = -1)+
  theme_void()

####

PAMPer.imp<-FindVariableFeatures(PAMPer.imp,selection.method = "vst",mean.cutoff = c(0.8,8),nfeatures = 600)
top10 <- rownames(PAMPer.imp)[grep("CER",rownames(PAMPer.imp))]
plot1 <- VariableFeaturePlot(PAMPer.imp)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) 
plot2
DotPlot(PAMPer.imp,features ="CER(20:1)",cols = c("white","red"))+RotatedAxis()
DE.lipdi.72h$class<-lipid_class$Class[match(rownames(DE.lipdi.72h),lipid_class$ï..Name)]
DE_Lipid_PE<-rownames(DE.lipdi.72h)[which(DE.lipdi.72h$p_val_adj<0.05 & DE.lipdi.72h$class=="PE")]
DE_Lipid_PC<-rownames(DE.lipdi.72h)[which(DE.lipdi.72h$p_val_adj<0.05 & DE.lipdi.72h$class=="PC")]
DE_Lipid_TAG<-rownames(DE.lipdi.72h)[which(DE.lipdi.72h$p_val_adj<0.05 & DE.lipdi.72h$class=="TAG")]
DE_Lipid_LPC<-rownames(DE.lipdi.72h)[which(DE.lipdi.72h$p_val_adj<0.05 & DE.lipdi.72h$class=="LPC")]
DE_Lipid_CE<-rownames(DE.lipdi.72h)[which(DE.lipdi.72h$p_val_adj<0.05 & DE.lipdi.72h$class=="CE")]

Major_class<-lipid_class$ï..Name[which(lipid_class$Class %in% c("TAG","DAG","PE","PC"))]
Major_class<-as.character(Major_class)
Minor_class<-setdiff(as.character(lipid_class$ï..Name),Major_class)
PE_class<-lipid_class$ï..Name[which(lipid_class$Class %in% c("PE"))]
PC_class<-lipid_class$ï..Name[which(lipid_class$Class %in% c("PC"))]
TAG_class<-lipid_class$ï..Name[which(lipid_class$Class %in% c("TAG"))]
DAG_class<-lipid_class$ï..Name[which(lipid_class$Class %in% c("DAG"))]

Choose_PE<-union(VariableFeatures(PAMPer.imp)[VariableFeatures(PAMPer.imp) %in% PE_class][1:44],DE_Lipid_PE)
Choose_PC<-union(VariableFeatures(PAMPer.imp)[VariableFeatures(PAMPer.imp) %in% PC_class][1:49],DE_Lipid_PC)
Choose_TAG<-union(VariableFeatures(PAMPer.imp)[VariableFeatures(PAMPer.imp) %in% TAG_class][1:63],DE_Lipid_TAG)
Choose_DAG<-VariableFeatures(PAMPer.imp)[VariableFeatures(PAMPer.imp) %in% DAG_class][1:40]

Choose_Lipid<-c(Choose_TAG,Choose_DAG,Choose_PC,Choose_PE,Minor_class)

Lipid.filter.data<-FetchData(PAMPer.imp,slot = "data",vars = Choose_Lipid)
cor_matrix<-cor(Lipid.filter.data)
cor_matrix[1:5,1:5]

cor_matrix[is.na(cor_matrix)]<-0
table(is.na(cor_matrix))
myPalette<-hue_pal()(length(unique(lipid_class$Class)))
names(myPalette)<-levels(lipid_class_stat$Class)
col_fun = colorRamp2(c(0, 0.5,1), c("slateblue2", "white", "firebrick2"))
ha = HeatmapAnnotation(Species = lipid_class$Class,col = list(Species=myPalette))
ht<-Heatmap(cor_matrix,name = "r",col =col_fun, show_row_names = T,
            show_column_names = F#,bottom_annotation = ha
            )
ht
write.csv(cor_matrix,file = "test.csv")
lipid_class2<-lipid_class[match(rownames(cor_matrix),lipid_class$ï..Name),]
write.csv(lipid_class2,file = "lipid_class2.csv")
lipid_class3<-lipid_class2
lipid_class3$Logfc<-DE.lipdi.72h$avg_logFC[match(lipid_class3$ï..Name,rownames(DE.lipdi.72h))]
lipid_class3$adjustp<-DE.lipdi.72h$p_val_adj[match(lipid_class3$ï..Name,rownames(DE.lipdi.72h))]
lipid_class3$Logfc[is.na(lipid_class3$Logfc)]<- 0
lipid_class3$adjustp[is.na(lipid_class3$adjustp)]<-1
lipid_class3$Lable<-"Unchanged"
lipid_class3$Lable[which(lipid_class3$adjustp<0.05 & lipid_class3$Logfc>0)]<-"Upregulated"
lipid_class3$Lable[which(lipid_class3$adjustp<0.05 & lipid_class3$Logfc<0)]<-"Downregulated"

write.csv(lipid_class3,file = "lipid_class3.csv")
####
saveRDS(lipid_class,file = "rdata/lipid_class.rdata")

DE.lipdi.72h$p_val[match(Lipid.sig,rownames(DE.lipdi.72h))]
