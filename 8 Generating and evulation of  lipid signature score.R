library(Seurat)
library(dplyr)
library(tidyr)
library(foreign)
library(stringr)
library(ggplot2)
library(EnhancedVolcano)

###Fig6A sFig6ABC
####
Lipid.sig<-covid2.lipid.pdata$PAMPer.ID[grep("PE",covid2.lipid.pdata$PAMPer.ID)]
Lipid.sig4<-Lipid.sig[c(1,2,6,4)]

###Preparation
Idents(PAMPer.imp)<-PAMPer.imp$TIME.POINT
PAMPer.imp.72h<-subset(PAMPer.imp,idents="72")
PAMPer.imp.72h<-ScaleData(PAMPer.imp.72h,features = rownames(PAMPer.imp.72h))
PAMPer.imp.72h$Lipid.sig.score<-apply(PAMPer.imp.72h@assays$RNA@scale.data[Lipid.sig,],2,mean)
PAMPer.imp.72h$Lipid.sig.score2<-apply(PAMPer.imp.72h@assays$RNA@scale.data[Lipid.sig2,],2,mean)
PAMPer.imp.72h$Lipid.sig.score3<-PAMPer.imp.72h@assays$RNA@scale.data[Lipid.sig3,]
PAMPer.imp.72h$Lipid.sig.score4<-apply(PAMPer.imp.72h@assays$RNA@scale.data[Lipid.sig4,],2,mean)

order<-names(sort(PAMPer.imp.72h$Lipid.sig.score))
PAMPer.imp.72h$clusters<-c(rep("Low",48),rep("Med",47),rep("High",47))[match(colnames(PAMPer.imp.72h),order)]
order<-names(sort(PAMPer.imp.72h$Lipid.sig.score4))
PAMPer.imp.72h$clusters2<-c(rep("Low",48),rep("Med",47),rep("High",47))[match(colnames(PAMPer.imp.72h),order)]
table(PAMPer.imp.72h$clusters,PAMPer.imp.72h$clusters2)


###volcano plot
PAMPer.imp.72h <- SetAssayData(object = PAMPer.imp.72h, slot = "data", new.data = as.matrix(log2(GetAssayData(object = PAMPer.imp.72h, slot = "counts")+1)))

Idents(PAMPer.imp.72h)<-PAMPer.imp.72h$outcome
DE.lipdi.72h<-FindMarkers(PAMPer.imp.72h,ident.1 = "Non-resolving",ident.2 = "Resolving",test.use = "LR",logfc.threshold = 0,latent.vars = c("INJURY.SEVERITY.SCORE.1","TRAUMATIC.BRAIN.INJURY","AGE","TREATMENT"))
DE.lipdi.72h<-DE.lipdi.72h[which(DE.lipdi.72h$avg_log2FC>0),]
DE.lipdi.72h$FC<-2^DE.lipdi.72h$avg_log2FC
write.csv(DE.lipdi.72h,file = "rdata/DE.lipdi.72h.csv")
EnhancedVolcano(DE.lipdi.72h,
                lab = rownames(DE.lipdi.72h),
                #selectLab = c(rownames(DE.lipdi.72h)[1:30]),
                x = 'avg_logFC',
                y = 'p_val_adj',
                drawConnectors = F,
                xlim = c(-3, 3),
                ylim = c(0,5),
                pCutoff = 0.01,
                FCcutoff = 0.4,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1, transcriptLabSize = 3,transcriptPointSize = 1)+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)
######correlation matrix

DE.lipids<-rownames(DE.lipdi.72h)[which( abs(DE.lipdi.72h$avg_log2FC) > 0.4 & DE.lipdi.72h$p_val_adj<0.01)] 
Test<-FetchData(PAMPer.imp,slot = "data",vars = union(DE.lipids,c(Lipid.sig)))
cor_matrix<-cor(Test,method = "spearman")

col_fun = colorRamp2(c(-1,0,1), c("Blue","white", "firebrick2"))
ht<-Heatmap(cor_matrix,name = "r",col =col_fun,
            cluster_columns = T,cluster_rows = T,show_row_names = T,show_column_names = F,
            width = unit(15, "cm"),height = unit(15, "cm")
            
)
ht

###UMAP in global pattern
PAMPer.imp$clusters<-PAMPer.imp.72h$clusters[match(PAMPer.imp$`PAMPER ID NUMBER`,PAMPer.imp.72h$`PAMPER ID NUMBER`)]
PAMPer.imp$clusters2<-PAMPer.imp.72h$clusters2[match(PAMPer.imp$`PAMPER ID NUMBER`,PAMPer.imp.72h$`PAMPER ID NUMBER`)]

table(PAMPer.imp$clusters2)
DimPlot(PAMPer.imp, reduction = "umap", label = F,pt.size = 3,ncol = 4,group.by = "clusters2",split.by = "outcome_time")+ ggplot2::theme_bw()+ggplot2::theme(legend.position = "bottom")+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)
####
PAMPer.imp$Lipid.sig.score<-apply(PAMPer.imp@assays$RNA@scale.data[Lipid.sig,],2,mean)
PAMPer.imp$Lipid.sig.score4<-apply(PAMPer.imp@assays$RNA@scale.data[Lipid.sig4,],2,mean)

FeaturePlot(PAMPer.imp,features = "Lipid.sig.score2",cols = c("lightgrey","red"),pt.size = 2)+ ggplot2::theme_bw()+ggplot2::theme(legend.position = "bottom")+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)
myPalette3<-c("#F8766D","#00BA38","#619CFF","white")
DimPlot(PAMPer.imp, reduction = "umap", label = F,pt.size = 2,ncol = 4,group.by = "clusters")+ ggplot2::theme_bw()+ggplot2::theme(legend.position = "bottom")+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)+scale_color_manual(values = myPalette3)

FeaturePlot(PAMPer.imp.72h,features = "Lipid.sig.score2",cols = c("lightgrey","red"),pt.size = 2,max.cutoff = 1)+ ggplot2::theme_bw()+ggplot2::theme(legend.position = "bottom")+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)
myPalette3<-c("#F8766D","#00BA38","#619CFF","white")
DimPlot(PAMPer.imp.72h, reduction = "umap", label = F,pt.size = 2,ncol = 4,group.by = "clusters")+ ggplot2::theme_bw()+ggplot2::theme(legend.position = "bottom")+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)+scale_color_manual(values = myPalette3)

##

###Applied to other 3 datasets
MM.imp$Lipid.sig.score<-apply(MM.imp@assays$RNA@scale.data[Lipid.sig,],2,mean)
covid.imp$Lipid.sig.score<-apply(covid.imp@assays$RNA@scale.data[Lipid.sig,],2,mean)
covid2.imp$Lipid.sig.score<-apply(covid2.imp@assays$RNA@scale.data[Lipid.sig,],2,mean)

MM.imp$Lipid.sig.score4<-apply(MM.imp@assays$RNA@scale.data[Lipid.sig4,],2,mean)
covid.imp$Lipid.sig.score4<-apply(covid.imp@assays$RNA@scale.data[Lipid.sig4,],2,mean)
covid2.imp$Lipid.sig.score4<-apply(covid2.imp@assays$RNA@scale.data[Lipid.sig4,],2,mean)
