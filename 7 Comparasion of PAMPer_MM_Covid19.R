###
library(readxl)
library(Seurat)
library(ComplexHeatmap)
library(dplyr)
library(pROC)
library(MASS)
library(ComplexHeatmap)
library(circlize)
library(dendsort)
library(ggplot2)
library(scales)

###Fig5 sFig5
MM.lipid.pdata$PAMPer_Lipid<-as.character(MM.lipid.pdata$PAMPer_Lipid)
covid.lipid.pdata$PAMPer_Lipid<-as.character(covid.lipid.pdata$PAMPer_Lipid)
covid2.lipid.pdata$PAMPer.ID<-as.character(covid2.lipid.pdata$PAMPer.ID)

###heatmap for common metabolites of PAMPer
PAMPer.heatmap<-FetchData(PAMPer.imp,slot = "scale.data",vars = c(covid2.lipid.pdata$PAMPer.ID,"outcome_time"))
PAMPer.heatmap2<-PAMPer.heatmap %>% group_by(outcome_time) %>% summarise_each(funs(mean)) 
timepoint<-PAMPer.heatmap2$outcome_time
PAMPer.heatmap2<-PAMPer.heatmap2[,-1]
PAMPer.heatmap2<-t(PAMPer.heatmap2)
colnames(PAMPer.heatmap2)<-timepoint

myPalette1<-hue_pal()(length(unique(MM.lipid.pdata$SUB_PATHWAY)))
names(myPalette1)<-unique(MM.lipid.pdata$SUB_PATHWAY)

col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))
ha = rowAnnotation(Class= MM.lipid.pdata$SUB_PATHWAY[match(rownames(PAMPer.heatmap2),MM.lipid.pdata$BIOCHEMICAL)],
                   col = list(Class=myPalette1))
ht<-Heatmap(PAMPer.heatmap2,name = "Exp",col =col_fun,
            show_row_names = T,cluster_columns = F,cluster_rows = F,
            right_annotation = ha,row_split = MM.lipid.pdata$SUB_PATHWAY[match(rownames(PAMPer.heatmap2),MM.lipid.pdata$BIOCHEMICAL)],
            row_gap = unit(3, "mm"),width = unit(4, "cm"),height = unit(10, "cm"))
ht

###heatmap for Common metabolites of MM
MM.heatmap<-FetchData(MM.imp,slot = "scale.data",vars = c(covid2.lipid.pdata$PAMPer.ID,"GROUP.TIME"))
MM.heatmap2<-MM.heatmap %>% group_by(GROUP.TIME) %>% summarise_each(funs(mean)) 
timepoint<-MM.heatmap2$GROUP.TIME
MM.heatmap2<-MM.heatmap2[,-1]
MM.heatmap2<-t(MM.heatmap2)
colnames(MM.heatmap2)<-timepoint


myPalette1<-hue_pal()(length(unique(MM.lipid.pdata$SUB_PATHWAY)))
names(myPalette1)<-unique(MM.lipid.pdata$SUB_PATHWAY)


col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))
ha = rowAnnotation(Class= MM.lipid.pdata$SUB_PATHWAY[match(rownames(MM.heatmap2),MM.lipid.pdata$BIOCHEMICAL)],
                   col = list(Class=myPalette1))
ht<-Heatmap(MM.heatmap2,name = "Exp",col =col_fun,
            show_row_names = T,cluster_columns = F,cluster_rows = F,
            right_annotation = ha,row_split = MM.lipid.pdata$SUB_PATHWAY[match(rownames(MM.heatmap2),MM.lipid.pdata$BIOCHEMICAL)],
            row_gap = unit(3, "mm"),width = unit(4, "cm"),height = unit(10, "cm"))
ht

###
###heatmap for Common metabolites of COVID.1
covid.heatmap<-FetchData(covid.imp,slot = "scale.data",vars = c(covid2.lipid.pdata$PAMPer.ID,"GROUP.TIME"))
covid.heatmap2<-covid.heatmap %>% group_by(GROUP.TIME) %>% summarise_each(funs(mean)) 
timepoint<-covid.heatmap2$GROUP.TIME
covid.heatmap2<-covid.heatmap2[,-1]
covid.heatmap2<-t(covid.heatmap2)
colnames(covid.heatmap2)<-timepoint

myPalette1<-hue_pal()(length(unique(MM.lipid.pdata$SUB_PATHWAY)))
names(myPalette1)<-unique(MM.lipid.pdata$SUB_PATHWAY)


col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))
ha = rowAnnotation(Class= MM.lipid.pdata$SUB_PATHWAY[match(rownames(covid.heatmap2),MM.lipid.pdata$BIOCHEMICAL)],
                   col = list(Class=myPalette1))
ht<-Heatmap(covid.heatmap2,name = "Exp",col =col_fun,
            show_row_names = T,cluster_columns = F,cluster_rows = F,
            right_annotation = ha,row_split = MM.lipid.pdata$SUB_PATHWAY[match(rownames(covid.heatmap2),MM.lipid.pdata$BIOCHEMICAL)]
            ,row_gap = unit(3, "mm"),width = unit(4, "cm"),height = unit(10, "cm")
)
ht
####heatmap for Common metabolites of COVID.2

covid2.heatmap<-FetchData(covid2.imp,slot = "scale.data",vars = c(covid2.lipid.pdata$PAMPer.ID,"GROUP.TIME"))
covid2.heatmap2<-covid2.heatmap %>% group_by(GROUP.TIME) %>% summarise_each(funs(mean)) 
timepoint<-covid2.heatmap2$GROUP.TIME
covid2.heatmap2<-covid2.heatmap2[,-1]
covid2.heatmap2<-t(covid2.heatmap2)
colnames(covid2.heatmap2)<-timepoint

myPalette1<-hue_pal()(length(unique(MM.lipid.pdata$SUB_PATHWAY)))
names(myPalette1)<-unique(MM.lipid.pdata$SUB_PATHWAY)


col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))
ha = rowAnnotation(Class= MM.lipid.pdata$SUB_PATHWAY[match(rownames(covid2.heatmap2),MM.lipid.pdata$BIOCHEMICAL)],
                   col = list(Class=myPalette1))
ht<-Heatmap(covid2.heatmap2,name = "Exp",col =col_fun,
            show_row_names = T,cluster_columns = F,cluster_rows = F,
            right_annotation = ha,row_split = MM.lipid.pdata$SUB_PATHWAY[match(rownames(covid2.heatmap2),MM.lipid.pdata$BIOCHEMICAL)]
            ,row_gap = unit(3, "mm"),width = unit(4, "cm"),height = unit(10, "cm")
)
ht

#####Statistics for all common metabolites in 4 heatmap
PAMPer.imp.common.meta<-subset(PAMPer.imp,features  = covid2.lipid.pdata$PAMPer.ID)
table(PAMPer.imp.common.meta$outcome_time)
Idents(PAMPer.imp.common.meta)<-PAMPer.imp.common.meta$outcome_time
PAMPer.common.meta.stat<-FindMarkers(PAMPer.imp.common.meta,ident.1 = "Non-resolving_72",ident.2 = "Resolving_72",logfc.threshold = 0)
PAMPer.common.meta.stat<-PAMPer.common.meta.stat[covid2.lipid.pdata$PAMPer.ID,]
write.csv(PAMPer.common.meta.stat[,c(1:2,5)],file = "Rdata/PAMPer.common.meta.stat.csv")

MM.imp.common.meta<-subset(MM.imp,features  = covid2.lipid.pdata$PAMPer.ID)
table(MM.imp.common.meta$GROUP.TIME)
Idents(MM.imp.common.meta)<-MM.imp.common.meta$GROUP.TIME
MM.common.meta.stat<-FindMarkers(MM.imp.common.meta,ident.1 = "Non-resolving_D2-D5",ident.2 = "Resolving_D2-D5",logfc.threshold = 0)
MM.common.meta.stat<-MM.common.meta.stat[covid2.lipid.pdata$PAMPer.ID,]
write.csv(MM.common.meta.stat[,c(1:2,5)],file = "Rdata/MM.common.meta.stat.csv")

covid.imp.common.meta<-subset(covid.imp,features  = covid2.lipid.pdata$PAMPer.ID)
table(covid.imp.common.meta$outcome)
Idents(covid.imp.common.meta)<-covid.imp.common.meta$outcome
covid.common.meta.stat<-FindMarkers(covid.imp.common.meta,ident.1 = "Severe Covid",ident.2 = "Mild-Mod Covid",logfc.threshold = 0)
covid.common.meta.stat<-covid.common.meta.stat[covid2.lipid.pdata$PAMPer.ID,]
write.csv(covid.common.meta.stat[,c(1:2,5)],file = "Rdata/covid.common.meta.stat.csv")

covid2.imp.common.meta<-subset(covid2.imp,features  = covid2.lipid.pdata$PAMPer.ID)
table(covid2.imp.common.meta$GROUP.TIME)
Idents(covid2.imp.common.meta)<-covid2.imp.common.meta$GROUP.TIME
covid2.common.meta.stat<-FindMarkers(covid2.imp.common.meta,ident.1 = c("SEVERE","MODERATE"),ident.2 = c("MILD","CONTROL"),logfc.threshold = 0)
covid2.common.meta.stat<-covid2.common.meta.stat[covid2.lipid.pdata$PAMPer.ID,]
write.csv(covid2.common.meta.stat[,c(1:2,5)],file = "Rdata/covid2.common.meta.stat.csv")
####Statistics for all common metabolites in 4 heatmap
PAMPer.imp.common.meta<-subset(PAMPer.imp,features  = as.character(MM.lipid.pdata$PAMPer_Lipid))
table(PAMPer.imp.common.meta$outcome_time)
Idents(PAMPer.imp.common.meta)<-PAMPer.imp.common.meta$outcome_time
PAMPer.common.meta.stat.75<-FindMarkers(PAMPer.imp.common.meta,ident.1 = "Non-resolving_72",ident.2 = "Resolving_72",logfc.threshold = 0)
write.csv(PAMPer.common.meta.stat.75[,c(1:2,5)],file = "Rdata/PAMPer.common.meta.stat.75.csv")

MM.imp.common.meta<-subset(MM.imp,features  = MM.lipid.pdata$PAMPer_Lipid)
table(MM.imp.common.meta$GROUP.TIME)
Idents(MM.imp.common.meta)<-MM.imp.common.meta$GROUP.TIME
MM.common.meta.stat.75<-FindMarkers(MM.imp.common.meta,ident.1 = "Non-resolving_D2-D5",ident.2 = "Resolving_D2-D5",logfc.threshold = 0)
write.csv(MM.common.meta.stat.75[,c(1:2,5)],file = "Rdata/MM.common.meta.stat.75.csv")
