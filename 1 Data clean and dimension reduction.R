library(dplyr)
library(tidyr)
library(gdata)
library(forcats)
library(foreign)
library(Seurat)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(dendsort)
library(FSA)
library(stringr)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

#Fig1CDE Fig2ABC sFig2 AB
###Data clean
####pdata in Total PAMPer trial
Pamper.pdata<-read.spss(file = "rawdata/ITT PAMPer data Josh variable sav.sav", use.value.label=TRUE, to.data.frame=TRUE)
Pamper.pdata<-Pamper.pdata %>% droplevels()
####Id transformation 
Pamper.ID<-read.csv(file = "rawdata/PAMP CRISMA_sperry.csv")
Pamper.pdata$PAMPID<-Pamper.ID$PAMPERID[match(Pamper.pdata$stnum,Pamper.ID$CRISMA.ID)]

Pamper.pdata$outcome<-"Non-resolving"
Pamper.pdata$outcome[which(Pamper.pdata$icu_los <7 & Pamper.pdata$t_30d_mort ==30)]<-"Resolving"
Pamper.pdata$outcome[which(Pamper.pdata$t_30d_mort <4)]<-"Early-Nonsurvivors"

Pamper.pdata$trans_icu_los<-Pamper.pdata$icu_los
Pamper.pdata$trans_icu_los[which(Pamper.pdata$outcome=="Early-Nonsurvivors")]<- -1
Pamper.pdata$trans_icu_los[which(Pamper.pdata$outcome=="Non-resolving" & Pamper.pdata$t_30d_mort<30)]<- 30
Pamper.pdata$trans_icu_los[which(Pamper.pdata$trans_icu_los>30)]<- 30

Pamper.pdata$trans_icu_los_censor<-1
Pamper.pdata$trans_icu_los_censor[which(Pamper.pdata$outcome=="Non-resolving" & Pamper.pdata$t_30d_mort<30)]<-0
Pamper.pdata$trans_icu_los_censor[which(Pamper.pdata$trans_icu_los>30)]<-0

table(Pamper.pdata$trans_icu_los,Pamper.pdata$trans_icu_los_censor)
Pamper.pdata$outcome[which(Pamper.pdata$t_30d_mort <4)]<-"Early-Nonsurvivors"
table(Pamper.pdata$outcome)

####pdata only in lipidomics 
pdata<-read.csv(file = "rawdata/pdata.csv",row.names = 1) %>%  t() %>% as.data.frame() 
pdata<-pdata %>% mutate(ID=rownames(pdata)) %>% filter(TREATMENT !="") %>% drop.levels()
level_key <- c('N/A' = "Baseline")
for (i in c("TREATMENT","TRAUMATIC BRAIN INJURY","ALIVE AT 30 DAYS","TIME POINT")) {
  pdata[,i]<-pdata[,i] %>% recode_factor( !!!level_key)
}
pdata$`PAMPER ID NUMBER`<-as.character(pdata$`PAMPER ID NUMBER`)
pdata$`PAMPER ID NUMBER`[which(pdata$TREATMENT=="Baseline")]<-as.character(pdata$`CLIENT SAMPLE ID`[which(pdata$TREATMENT=="Baseline")])
pdata$outcome<-Pamper.pdata$outcome[match(pdata$`PAMPER ID NUMBER`,Pamper.pdata$PAMPID)]
pdata$outcome[which(pdata$TREATMENT=="Baseline")]<-"Baseline"
pdata$outcome<-factor(pdata$outcome,levels = c("Baseline","Resolving","Non-resolving","Early-Nonsurvivors"))
pdata$AGEGROUP<-pdata$AGE %>% as.numeric.factor() %>% cut(breaks=c(0,40,60,200),labels = c("young","middle","old")) 
pdata$Severity<-pdata$`INJURY SEVERITY SCORE 1` %>% as.numeric.factor() %>% cut(breaks=c(0,14,24,75),labels = c("Mild","Mod","Severe")) %>% as.character()
pdata$Severity[is.na(pdata$Severity)]<-"Baseline"
pdata$Severity<-factor(pdata$Severity,levels = c("Baseline","Mild","Mod","Severe"))
pdata$TREATMENT<-pdata$TREATMENT %>% fct_relevel("Baseline","standard care","prehospital plasma")
pdata$`ALIVE AT 30 DAYS`<-pdata$`ALIVE AT 30 DAYS` %>% fct_relevel("Baseline","yes","no")

for (i in c("TREATMENT","TRAUMATIC BRAIN INJURY","ALIVE AT 30 DAYS","TIME POINT","Severity","AGEGROUP","outcome","PAMPER ID NUMBER")) {
  pdata[,paste(i,"time",sep = "_")]=interaction(pdata[,i],pdata[,"TIME POINT"],sep = "_",lex.order = F,drop = T) 
  
}
Filter_ID<-names(which(table(pdata[which(pdata$outcome %in% c("Non-resolving")),"PAMPER ID NUMBER"])<3))
pdata<-pdata %>% filter(!outcome_time %in% c("Early-Nonsurvivors_24","Early-Nonsurvivors_72")) %>% filter(!duplicated(`PAMPER ID NUMBER_time`)) %>% filter(`PAMPER ID NUMBER`!=Filter_ID) 
pdata$outcome_time<-pdata$outcome_time %>% droplevels()
rownames(pdata)<-pdata$ID

######Dimension reduction
Impdata<-read.csv(file = "rawdata/pdata.imp.csv",header = T,row.names = 1)
Impdata[1:5,1:5]
Impdata<-Impdata[,rownames(pdata)]
PAMPer.imp<-CreateSeuratObject(counts = Impdata,meta.data = pdata)
PAMPer.imp <- FindVariableFeatures(PAMPer.imp, selection.method = "vst", nfeatures = nrow(PAMPer.imp))
PAMPer.imp <- ScaleData(PAMPer.imp,do.scale = T,do.center = T)
PAMPer.imp <- RunPCA(PAMPer.imp, npcs = 30, verbose = FALSE)
PAMPer.imp <- RunUMAP(PAMPer.imp, reduction = "pca", dims = 1:20)
myPalette<-hue_pal()(length(unique(PAMPer.imp$`TIME POINT`)))
myPalette2<-hue_pal()(length(unique(PAMPer.imp$outcome_time)))
myPalette2<-c("#F8766D","green","royalblue","yellow2","green2","skyblue2","green4","royalblue4")
myPalette3<-c("#F8766D","#00BA38","#619CFF","white")
PAMPer.imp$TIME.POINT<-PAMPer.imp$`TIME POINT`
DimPlot(PAMPer.imp, reduction = "umap", label = T,pt.size = 2,ncol = 2,group.by = "TIME.POINT")+ ggplot2::theme_bw()+ggplot2::theme(legend.position = "bottom")+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)+scale_color_manual(values = myPalette)
DimPlot(PAMPer.imp, reduction = "umap", label = T,pt.size = 2,ncol = 2,group.by = "outcome_time")+ ggplot2::theme_bw()+ggplot2::theme(legend.position = "bottom")+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)+scale_color_manual(values = myPalette2)
DimPlot(PAMPer.imp, reduction = "umap", label = F,pt.size = 0.5,ncol = 4,group.by = "outcome_time",split.by = "outcome_time")+ ggplot2::theme_bw()+ggplot2::theme(legend.position = "bottom")+theme(panel.grid.major=element_blank(),panel.grid.minor =element_blank())+theme(aspect.ratio=1)
####Total.conc & boxplot
PAMPer.imp$Total.conc<-apply(PAMPer.imp@assays$RNA@data,2,sum)
myPalette<-hue_pal()(length(unique(PAMPer.imp$outcome)))

ggplot(PAMPer.imp@meta.data, aes(x=TIME.POINT, y=Total.conc,fill=TIME.POINT))+  geom_boxplot()+theme_classic()+theme(aspect.ratio=1,legend.position = "bottom")+theme(axis.text.x = element_text(angle = 0))+ylim(0,3000)+scale_color_manual(values = myPalette)
###Normality check
shapiro_test( PAMPer.imp$Total.conc[which(PAMPer.imp$`TIME POINT`== "Baseline")])
####kruskal.test
kruskal.test(Total.conc ~ TIME.POINT,data = PAMPer.imp@meta.data)
### Dunn test
dunnTest(Total.conc ~ TIME.POINT,
              data=PAMPer.imp@meta.data,
              method="bh")    # Can adjust p-values;
###Group based heatmap

PAMPer.heatmap<-FetchData(PAMPer.imp,slot = "scale.data",vars = c(VariableFeatures(PAMPer.imp)[1:900],"outcome_time"))
PAMPer.heatmap2<-PAMPer.heatmap %>% group_by(outcome_time) %>% summarise_each(funs(mean)) 
timepoint<-PAMPer.heatmap2$outcome_time
PAMPer.heatmap2<-PAMPer.heatmap2[,-1]
PAMPer.heatmap2<-t(PAMPer.heatmap2)
colnames(PAMPer.heatmap2)<-timepoint


col_fun = colorRamp2(c(-1,0,1), c("Blue","white", "firebrick2"))

ht<-Heatmap(PAMPer.heatmap2,name = "Exp",col =col_fun,
            show_row_names = F,cluster_columns = F,width = unit(10, "cm"),height = unit(10, "cm"))
ht
##Sub-Group based heatmap
Idents(PAMPer.imp)<-PAMPer.imp$`TRAUMATIC BRAIN INJURY`
PAMPer.heatmap<-FetchData(PAMPer.imp %>% subset(idents=c("Baseline","yes")),slot = "scale.data",vars = c(rownames(PAMPer.imp),"outcome_time"))
PAMPer.heatmap2<-PAMPer.heatmap %>% group_by(outcome_time) %>% summarise_each(funs(mean)) 
timepoint<-PAMPer.heatmap2$outcome_time
PAMPer.heatmap2<-PAMPer.heatmap2[,-1]
PAMPer.heatmap2<-t(PAMPer.heatmap2)
colnames(PAMPer.heatmap2)<-timepoint


col_fun = colorRamp2(c(-1,0,1), c("Blue","white", "firebrick2"))

ht<-Heatmap(PAMPer.heatmap2,name = "Exp",col =col_fun,
            show_row_names = F,cluster_columns = F,width = unit(10, "cm"),height = unit(10, "cm"))
ht
###Preparing for HC data
Idents(PAMPer.imp)<-PAMPer.imp$outcome
table(PAMPer.imp$outcome)
HC.data<-t(FetchData(subset(PAMPer.imp,idents = c("Baseline")),vars = rownames(PAMPer.imp) ))
rownames(HC.data)[1:10]
rownames(HC.data)[grep("TAG",rownames(HC.data))]<-rownames(HC.data)[grep("TAG",rownames(HC.data))] %>% str_sub(1,7) %>% str_replace("TAG","TAG ")
#rownames(HC.data)<-rownames(HC.data) %>% str_replace("A","")
rownames(HC.data)[-c(grep("^PE",rownames(HC.data)),grep("^PC",rownames(HC.data)),grep("DAG",rownames(HC.data)),grep("PI",rownames(HC.data)))]<-rownames(HC.data)[-c(grep("^PE",rownames(HC.data)),grep("^PC",rownames(HC.data)),grep("DAG",rownames(HC.data)),grep("PI",rownames(HC.data)))] %>% str_replace("[(]"," ") %>% str_replace("[)]","")

SM.name<-data.frame(raw=  rownames(HC.data)[grep("SM",rownames(HC.data))],
                     B= rownames(HC.data)[grep("SM",rownames(HC.data))] %>% str_sub(1,2),
                    C= rownames(HC.data)[grep("SM",rownames(HC.data))] %>% str_sub(4,5),
                    D= rownames(HC.data)[grep("SM",rownames(HC.data))] %>% str_sub(7))
SM.name$B<-SM.name$B %>% str_replace("SM","SM d")
SM.name$C<-SM.name$C %>% as.numeric.factor() %>% +18
SM.name$D<-SM.name$D %>% as.numeric.factor() %>% +1
rownames(HC.data)[grep("SM",rownames(HC.data))]<-SM.name$B %>% paste(SM.name$C,sep = "") %>% paste(SM.name$D,sep = ":")

CER.name<-data.frame(raw=  rownames(HC.data)[grep("^CER",rownames(HC.data))],
                    B= rownames(HC.data)[grep("^CER",rownames(HC.data))] %>% str_sub(1,3),
                    C= rownames(HC.data)[grep("^CER",rownames(HC.data))] %>% str_sub(5,6),
                    D= rownames(HC.data)[grep("^CER",rownames(HC.data))] %>% str_sub(8))
CER.name$B<-CER.name$B %>% str_replace("CER","CER d")
CER.name$C<-CER.name$C %>% as.numeric.factor() %>% +18
CER.name$D<-CER.name$D %>% as.numeric.factor() %>% +1
rownames(HC.data)[grep("^CER",rownames(HC.data))]<-CER.name$B %>% paste(CER.name$C,sep = "") %>% paste(CER.name$D,sep = ":")

HCER.name<-data.frame(raw=  rownames(HC.data)[grep("^HCER",rownames(HC.data))],
                     B= rownames(HC.data)[grep("^HCER",rownames(HC.data))] %>% str_sub(1,4),
                     C= rownames(HC.data)[grep("^HCER",rownames(HC.data))] %>% str_sub(6,7),
                     D= rownames(HC.data)[grep("^HCER",rownames(HC.data))] %>% str_sub(9))
HCER.name$B<-HCER.name$B %>% str_replace("HCER","HexCer d")
HCER.name$C<-HCER.name$C %>% as.numeric.factor() %>% +18
HCER.name$D<-HCER.name$D %>% as.numeric.factor() %>% +1
rownames(HC.data)[grep("^HCER",rownames(HC.data))]<-HCER.name$B %>% paste(HCER.name$C,sep = "") %>% paste(HCER.name$D,sep = ":")

DCER.name<-data.frame(raw=  rownames(HC.data)[grep("^DCER",rownames(HC.data))],
                      B= rownames(HC.data)[grep("^DCER",rownames(HC.data))] %>% str_sub(1,4),
                      C= rownames(HC.data)[grep("^DCER",rownames(HC.data))] %>% str_sub(6,7),
                      D= rownames(HC.data)[grep("^DCER",rownames(HC.data))] %>% str_sub(9))
DCER.name$B<-DCER.name$B %>% str_replace("DCER","CerOH d")
DCER.name$C<-DCER.name$C %>% as.numeric.factor() %>% +18
DCER.name$D<-DCER.name$D %>% as.numeric.factor() %>% +1
rownames(HC.data)[grep("^DCER",rownames(HC.data))]<-DCER.name$B %>% paste(DCER.name$C,sep = "") %>% paste(DCER.name$D,sep = ":")

rownames(HC.data)[grep("^CerOH",rownames(HC.data))]

write.csv(HC.data,file = "HC.data.csv")

