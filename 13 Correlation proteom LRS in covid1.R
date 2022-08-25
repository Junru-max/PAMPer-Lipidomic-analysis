####
library(readxl)
library(Seurat)
library(ComplexHeatmap)
library(dplyr)
library(pROC)
library(MASS)
library(circlize)
library(dendsort)
library(ggplot2)
library(scales)
library(stringr)
library("clusterProfiler")
library("DOSE")
library("org.Hs.eg.db")
library("enrichplot")
library("ReactomePA")
######Fig.7B sFig.7CDEF
###For lipid data
Idents(covid.imp)<-covid.imp$TIME
covid.imp.pa<-subset(covid.imp,idents=c("Before","<48h","D2-D5","D6_D14"))
covid.imp.pa<-ScaleData(covid.imp.pa,features = rownames(covid.imp.pa))
covid.imp.pa$Lipid.sig.score<-apply(covid.imp.pa@assays$RNA@scale.data[Lipid.sig,],2,mean)
covid.imp.pa$Lipid.sig.score2<-apply(covid.imp.pa@assays$RNA@scale.data[Lipid.sig2,],2,mean)

hist(covid.imp.pa$Lipid.sig.score2)
###For protein data
covid.pro<-read.csv(file = "rawdata/protein-omics_covid.csv",row.names = 1)
covid.pro.id<-data.frame(rawid=rownames(covid.pro),geneid=covid.pro$Gene.Symbol)
covid.pro<-covid.pro[,-1]
table(colnames(covid.pro) %in% covid.pdata$MS.ID.b)
covid.pro[is.na(covid.pro)]<-0
covid.pro<-covid.pro[,which(colnames(covid.pro) %in% covid.imp.pa$MS.ID.b)]

dim(covid.pro)
covid.pro<-as.data.frame(scale(t(covid.pro)))

covid.pro$Lipid.sig.score<-covid.imp.pa$Lipid.sig.score2[match(rownames(covid.pro),covid.imp.pa$MS.ID.b)]


cor_matrix<-cor(covid.pro,method = "spearman")
hist(cor_matrix[,"Lipid.sig.score"],breaks = 100)
cor_matrix[as.character(covid.pro.id$rawid[grep("APO",covid.pro.id$geneid)]),"Lipid.sig.score"]

Pos.pro<-cor_matrix[,"Lipid.sig.score"][which(cor_matrix[,"Lipid.sig.score"]> 0.3)]
Neg.pro<-cor_matrix[,"Lipid.sig.score"][which(cor_matrix[,"Lipid.sig.score"]< -0.3)]

###
Pos.pro.id<-covid.pro.id$geneid[which(covid.pro.id$rawid %in% names(Pos.pro))]
Neg.pro.id<-covid.pro.id$geneid[which(covid.pro.id$rawid %in% names(Neg.pro))]
Neg.pro.id<-Neg.pro.id[!is.na(Neg.pro.id)]

Pos.pro.id<-as.character(Pos.pro.id)
Neg.pro.id<-as.character(Neg.pro.id)

names(Pos.pro.id)<-as.character(Pos.pro.id)
names(Neg.pro.id)<-as.character(Neg.pro.id)

Pos.pro.id1<-Pos.pro.id[grep(",",Pos.pro.id)]
Pos.id1.name<-names(Pos.pro.id1)
Pos.pro.id2<-Pos.pro.id[setdiff(1:length(Pos.pro.id1),match(Pos.pro.id1,Pos.pro.id))]

Neg.pro.id1<-Neg.pro.id[grep(",",Neg.pro.id)]
Neg.id1.name<-names(Neg.pro.id1)
Neg.pro.id2<-Neg.pro.id[setdiff(1:length(Neg.pro.id1),match(Neg.pro.id1,Neg.pro.id))]



Pos.pro.id1<-Pos.pro.id1 %>% str_extract(".+?\\,") %>% str_sub(1,-2)
names(Pos.pro.id1)<-Pos.id1.name
Pos.pro.id<-c(Pos.pro.id1,Pos.pro.id2)

Neg.pro.id1<-Neg.pro.id1 %>% str_extract(".+?\\,") %>% str_sub(1,-2)
names(Neg.pro.id1)<-Neg.id1.name
Neg.pro.id<-c(Neg.pro.id1,Neg.pro.id2)
########Filter pos pro and neg pro

Cor.pdata<-data.frame(Genesymbol_C= c(names(c(Pos.pro.id,Neg.pro.id))),Genesymbol_S=c(Pos.pro.id,Neg.pro.id))
Cor.pdata$proid<-covid.pro.id$rawid[match(Cor.pdata$Genesymbol_C,covid.pro.id$geneid)]

####Pos pathway
Hs<-org.Hs.eg.db
ID_e_s_all_0<-select(Hs,
                     keys = Pos.pro.id,
                     columns = c("ENTREZID","SYMBOL"),
                     keytype = "SYMBOL")
ID_e_s_all_0<-ID_e_s_all_0[!is.na(ID_e_s_all_0$ENTREZID),]
PA<-enrichPathway(as.character(ID_e_s_all_0$ENTREZID), organism = "human", pvalueCutoff = 0.05,
              pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 2,
              maxGSSize = 500, readable = TRUE)
excludeCCdown <- data.frame("ID" = PA@result$ID, "name" = PA@result$Description,"gene"= PA@result$geneID)
selected.pathway <- c("R-HSA-166658","R-HSA-114608","R-HSA-6798695","R-HSA-2173782","R-HSA-948021")
excludeCCdown <- PA$ID[! excludeCCdown$ID %in% selected.pathway]
PA.tmp2 <- dropGO(PA, term = excludeCCdown)
dotplot(PA.tmp2,showCategory=15)


#####Neg pathway
Hs<-org.Hs.eg.db
ID_e_s_all_0<-select(Hs,
                     keys = Neg.pro.id,
                     columns = c("ENTREZID","SYMBOL"),
                     keytype = "SYMBOL")
ID_e_s_all_0<-ID_e_s_all_0[!is.na(ID_e_s_all_0$ENTREZID),]
PA<-enrichPathway(as.character(ID_e_s_all_0$ENTREZID), organism = "human", pvalueCutoff = 0.05,
                  pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 2,
                  maxGSSize = 500, readable = TRUE)
PA.result<-data.frame(PA)

excludeCCdown <- data.frame("ID" = PA@result$ID, "name" = PA@result$Description,"gene"= PA@result$geneID)
selected.pathway <- c("R-HSA-381426","R-HSA-8957275","R-HSA-76005","R-HSA-2173782","R-HSA-8964058","R-HSA-2168880")
excludeCCdown <- PA$ID[! excludeCCdown$ID %in% selected.pathway]
PA.tmp2 <- dropGO(PA, term = excludeCCdown)
dotplot(PA.tmp2,showCategory=15)


###correlation heatmap neg
select.gene<-c(Neg.pro.id)
select.gene<-c(Pos.pro.id)


covid.imp.line.data<-covid.pro[,as.character(Cor.pdata$proid)[Cor.pdata$Genesymbol_S %in% select.gene]]
colnames(covid.imp.line.data)<-Cor.pdata$Genesymbol_S[match(colnames(covid.imp.line.data),Cor.pdata$proid)]
covid.imp.line.data$Lipid.sig.score<-covid.imp.pa$Lipid.sig.score2[match(rownames(covid.imp.line.data),covid.imp.pa$MS.ID.b)]

cor_matrix<-cor(covid.imp.line.data,method = "spearman")
cor_matrix[is.na(cor_matrix)]<-0
table(is.na(cor_matrix))
ht<-Heatmap(cor_matrix,name = "r",col =col_fun,
            cluster_columns = T,cluster_rows = T,show_row_names = T,show_column_names = F,
            width = unit(10, "cm"),height = unit(10, "cm"))


ht
###correlation heatmap pos
select.gene<-c(Pos.pro.id)

covid.imp.line.data<-covid.pro[,as.character(Cor.pdata$proid)[Cor.pdata$Genesymbol_S %in% select.gene]]
colnames(covid.imp.line.data)<-Cor.pdata$Genesymbol_S[match(colnames(covid.imp.line.data),Cor.pdata$proid)]
covid.imp.line.data$Lipid.sig.score<-covid.imp.pa$Lipid.sig.score2[match(rownames(covid.imp.line.data),covid.imp.pa$MS.ID.b)]

cor_matrix<-cor(covid.imp.line.data,method = "spearman")
cor_matrix[is.na(cor_matrix)]<-0
table(is.na(cor_matrix))

number_k=25
Label<-c(sample(1:nrow(cor_matrix),size = number_k),124)
ma = rowAnnotation(foo= anno_mark(at = Label,
                                  labels= rownames(cor_matrix)[Label]))

col_fun = colorRamp2(c(-1,0,1), c("Blue","white", "firebrick2"))
ht<-Heatmap(cor_matrix,name = "r",col =col_fun,
            cluster_columns = T,cluster_rows = T,show_row_names = F,show_column_names = F,
            width = unit(10, "cm"),height = unit(10, "cm")
            ,right_annotation = ma
)
ht
###



####Selected heatmap
Acute.phase.gene<-c("SERPINA1","ITIH4","SERPINA3","ORM1","ORM2","SAA4","CRP","HP","SAA1","LBP","SAA2")
Complement.gene<-c("C4BPA","CFI","SERPING1","CRP","C7","C6","C2","C8A","C9","C8G")
Immunoglobulin.gene<-c("IGLC6","IGHV1-3","IGHV3-23","IGKV2D-30","IGHV3-21","IGKV3D-20","IGLV1-47")
IGF.gene<-c("ALB","FN1","IGF2","IGFALS","IGFBP3","IGFBP5","ITIH2")

select.gene<-c(Acute.phase.gene,Complement.gene,Immunoglobulin.gene,IGF.gene)

covid.imp.line.data<-covid.pro[,as.character(Cor.pdata$proid)[Cor.pdata$Genesymbol_S %in% select.gene]]
colnames(covid.imp.line.data)<-Cor.pdata$Genesymbol_S[match(colnames(covid.imp.line.data),Cor.pdata$proid)]
covid.imp.line.data$Lipid.sig.score<-covid.imp.pa$Lipid.sig.score2[match(rownames(covid.imp.line.data),covid.imp.pa$MS.ID.b)]

cor_matrix<-cor(covid.imp.line.data,method = "spearman")

cor_matrix[is.na(cor_matrix)]<-0

table(is.na(cor_matrix))

col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))
anno_df = data.frame(Acute_phase = rep(0,length(colnames(cor_matrix))),
                     Complement = rep(0,length(colnames(cor_matrix))),
                     Immunoglobulin = rep(0,length(colnames(cor_matrix))),
                     IGF=rep(0,length(colnames(cor_matrix))))
rownames(anno_df)<-rownames(cor_matrix)
anno_df$Acute_phase[which(rownames(cor_matrix) %in% Acute.phase.gene)]<-1
anno_df$Complement[which(rownames(cor_matrix) %in% Complement.gene)]<-1
anno_df$Immunoglobulin[which(rownames(cor_matrix) %in% Immunoglobulin.gene)]<-1

anno_df$IGF[which(rownames(cor_matrix) %in% IGF.gene)]<-1


ha = rowAnnotation(df = anno_df,
                       col = list(Acute_phase = c("0" = "white", "1" = "red" ),
                                  Complement = c("0" = "white", "1" = "green"),
                                  Immunoglobulin = c("0" = "white", "1" = "Orange"),
                                  IGF=c("0"="white","1"="blue")
                       )
)

ht<-Heatmap(cor_matrix,name = "z-score",col =col_fun,
            cluster_columns = T,cluster_rows = T,show_row_names = T,show_column_names = F
            ,right_annotation = ha
            ,
            width = unit(10, "cm"),height = unit(10, "cm")
               
)


ht
