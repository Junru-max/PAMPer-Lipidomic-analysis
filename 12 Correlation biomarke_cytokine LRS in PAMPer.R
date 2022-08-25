library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(dendsort)
##Fig.7A sFig.7AB
######

###LIPID
Idents(PAMPer.imp.72h)<-PAMPer.imp.72h$TIME.POINT
PAMPer.imp.72h$PAID<-PAMPer.imp.72h$`PAMPER ID NUMBER`
PAMPer.lipid.72h<-FetchData(subset(PAMPer.imp.72h,ident=72),vars = c("PAID","Lipid.sig.score2"))
####
###Biomarkers
Idents(biomarker)<-biomarker$timepoint
PAMPer.biomarker.0h<-FetchData(subset(biomarker,ident="0"),vars = c("PAID",rownames(biomarker)))
PAMPer.biomarker.0h[,2:8]<-log2(PAMPer.biomarker.0h[,2:8]+1)
PAMPer.biomarker.24h<-FetchData(subset(biomarker,ident="24"),vars = c("PAID",rownames(biomarker)))
PAMPer.biomarker.24h[,2:8]<-log2(PAMPer.biomarker.24h[,2:8]+1)
###
###Biomarkers

Idents(Luminex)<-Luminex$timepoint
PAMPer.Luminex.0h<-FetchData(subset(Luminex,ident="0"),vars = c("PAID",rownames(Luminex)))
PAMPer.Luminex.0h[,2:22]<-log2(PAMPer.Luminex.0h[,2:22]+1)
PAMPer.Luminex.24h<-FetchData(subset(Luminex,ident="24"),vars = c("PAID",rownames(Luminex)))
PAMPer.Luminex.24h[,2:22]<-log2(PAMPer.Luminex.24h[,2:22]+1)
PAMPer.Luminex.72h<-FetchData(subset(Luminex,ident="72"),vars = c("PAID",rownames(Luminex)))
PAMPer.Luminex.72h[,2:22]<-log2(PAMPer.Luminex.72h[,2:22]+1)
###
PAMPer.luminex_lipid<-inner_join(PAMPer.biomarker.0h,PAMPer.Luminex.0h,by="PAID") %>% inner_join(PAMPer.lipid.72h,by="PAID")
PAMPer.luminex_lipid<-PAMPer.luminex_lipid[,-which(colnames(PAMPer.luminex_lipid) %in% c("PAID"))]
PAMPer.luminex_lipid<-as.data.frame(PAMPer.luminex_lipid)

##
### correlation heatmap 
cor_matrix<-cor(PAMPer.luminex_lipid,method = "spearman")
cor_matrix[1:5,1:5]
cor_matrix[is.na(cor_matrix)]<-0
cor_matrix[,"Lipid.sig.score2"]

### Prepare for annotation
Subset1<-c("IL.6","IL.8","IL.10","MCP.1","MIG","IP.10","TNFa")
Subset2<-c("IL.4","IL.5","IL.1b", "IL.2","IL.7", "IL.17A","GM.CSF")
Subset3<-c("IL.17E","IL.22","IL.9","IL.33","IL.21","IL.23.ng.ml.","IL.27")
EC<-c("sFLT.1.sVEGFR1-pg.mL","Syndecan.1-ng.ml","Thrombomodulin-ng.ml","suPAR-ng.mL","Cell.death-","S100A10-ng.mL")
Adipokine<-c("Adiponectin-ng.mL")

col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))
anno_df = data.frame(Subset1 = rep(0,length(colnames(cor_matrix))),
                     Subset2 = rep(0,length(colnames(cor_matrix))),
                     Subset3 = rep(0,length(colnames(cor_matrix))),
                     EC=rep(0,length(colnames(cor_matrix))),
                     Adipokine=rep(0,length(colnames(cor_matrix)))
                     )
rownames(anno_df)<-rownames(cor_matrix)
anno_df$Subset1[which(rownames(cor_matrix) %in% Subset1)]<-1
anno_df$Subset2[which(rownames(cor_matrix) %in% Subset2)]<-1
anno_df$Subset3[which(rownames(cor_matrix) %in% Subset3)]<-1
anno_df$EC[which(rownames(cor_matrix) %in% EC)]<-1
anno_df$Adipokine[which(rownames(cor_matrix) %in% Adipokine)]<-1


ha = rowAnnotation(df = anno_df,
                   col = list(Subset1 = c("0" = "white", "1" = hue_pal()(8)[1] ),
                              Subset2 = c("0" = "white", "1" = hue_pal()(8)[2]),
                              Subset3 = c("0" = "white", "1" = hue_pal()(8)[3]),
                              EC = c("0" = "white", "1" = hue_pal()(8)[5]),
                              Adipokine = c("0" = "white", "1" = hue_pal()(8)[6])
                   )
           )

ht<-Heatmap(cor_matrix,name = "r",col =col_fun,
            cluster_columns = T,cluster_rows = T,show_row_names = T,show_column_names = F
            ,right_annotation = ha,
            width = unit(15, "cm"),height = unit(15, "cm"),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.1f", cor_matrix[i, j]), x, y, gp = gpar(fontsize = 5))
            }
            
)
ht

#### 24
###
PAMPer.luminex_lipid<-inner_join(PAMPer.biomarker.24h,PAMPer.Luminex.24h,by="PAID") %>% inner_join(PAMPer.lipid.72h,by="PAID")
PAMPer.luminex_lipid<-PAMPer.luminex_lipid[,-which(colnames(PAMPer.luminex_lipid) %in% c("PAID"))]
PAMPer.luminex_lipid<-as.data.frame(PAMPer.luminex_lipid)

##
### correlation heatmap 
cor_matrix<-cor(PAMPer.luminex_lipid,method = "spearman")
cor_matrix[1:5,1:5]
cor_matrix[is.na(cor_matrix)]<-0

group<-rownames(cor_matrix)
group[which(group %in% Subset1)]<-"Subset1 Cytokines"
group[which(group %in% Subset2)]<-"Subset2 Cytokines"
group[which(group %in% Subset3)]<-"Subset3 Cytokines"
group[which(group %in% EC)]<-"EC injury biomarkers"
group[which(group %in% Adipokine)]<-"Adipokine"

group<-factor(group)
myPalette2<-hue_pal()(length(unique(group)))
names(myPalette2)<-levels(group)
ha = rowAnnotation(Class= group,
                   col = list(Class=myPalette2))

col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))

ht<-Heatmap(cor_matrix,name = "Exp",col =col_fun, show_row_names = T,cluster_columns = T,
            width = unit(15, "cm"),height = unit(15, "cm"),right_annotation = ha,show_row_dend = F,show_column_dend = T)
ht
###  72h
### correlation heatmap
PAMPer.lipid.72h$PAID
dim(PAMPer.Luminex.72h)
PAMPer.luminex_lipid<-PAMPer.Luminex.72h %>% inner_join(PAMPer.lipid.72h,by="PAID")
PAMPer.luminex_lipid<-PAMPer.luminex_lipid[,-which(colnames(PAMPer.luminex_lipid) %in% c("PAID"))]
PAMPer.luminex_lipid<-as.data.frame(PAMPer.luminex_lipid)

cor_matrix<-cor(PAMPer.luminex_lipid,method = "spearman")
cor_matrix[1:5,1:5]
cor_matrix[is.na(cor_matrix)]<-0

group<-rownames(cor_matrix)
group[which(group %in% Subset1)]<-"Subset1 Cytokines"
group[which(group %in% Subset2)]<-"Subset2 Cytokines"
group[which(group %in% Subset3)]<-"Subset3 Cytokines"

group<-factor(group)
myPalette2<-hue_pal()(6)[3:6]
names(myPalette2)<-levels(group)
ha = rowAnnotation(Class= group,
                   col = list(Class=myPalette2))

col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))

ht<-Heatmap(cor_matrix,name = "Exp",col =col_fun, show_row_names = T,cluster_columns = T,
            width = unit(15, "cm"),height = unit(15, "cm"),right_annotation = ha,show_row_dend = F,show_column_dend = T)
ht
