library(Seurat)
library(dplyr)
library(tidyr)
library(foreign)
library(stringr)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(dendsort)

####sFig 3BC
####Lipidomics_pdata
Idents(PAMPer.imp)<-PAMPer.imp$outcome
data_merge<-subset(PAMPer.imp,idents="Baseline",invert=T)
data_merge$PAID<-data_merge$`PAMPER ID NUMBER`
data_merge$PAID_TIME<-paste(data_merge$PAID,data_merge$TIME.POINT,sep = "_")

####Luminex_pdata
PAMPer.luminex<-read.csv(file = "rawdata/PAMPer luminex .csv",row.names = 1)
dim(PAMPer.luminex)
#PAMPer.luminex[is.na(PAMPer.luminex)]<-0
#PAMPer.luminex<-PAMPer.luminex[-which(apply(PAMPer.luminex, 1, sum)==0),]
timepoint<-rownames(PAMPer.luminex) %>% str_extract("(?<=\\-).+?(?=\\-)")
timepoint<- timepoint %>% gsub(pattern = "\\b0H\\b",replacement = "0HR") %>% gsub(pattern = "hr",replacement = "HR") %>% gsub(pattern = "HR",replacement = "") 
PAID<-rownames(PAMPer.luminex) %>% str_extract(".+?\\-") %>% str_extract(".+?\\-") %>% str_sub(5,-2)
Luminex.pdata<-data.frame(rownames =rownames(PAMPer.luminex),PAID=PAID, timepoint=timepoint,PAID_TIME=paste(PAID,timepoint,sep = "_"))
Luminex.pdata$outcome<-Pamper.pdata$outcome[match(Luminex.pdata$PAID,Pamper.pdata$PAMPID)]
Luminex.pdata$outcome_time<-paste(Luminex.pdata$outcome,Luminex.pdata$timepoint,sep = "_")
dim(Luminex.pdata)
Filter_ID<-names(which(table(as.character(Luminex.pdata[which(Luminex.pdata$outcome %in% c("Non-resolving","Resolving")),"PAID"]))<3))

Luminex.pdata<-Luminex.pdata %>% filter(!outcome_time %in% c("Early-Nonsurvivors_24","Early-Nonsurvivors_72"))  %>% filter(!PAID %in% Filter_ID) 

table(duplicated(Luminex.pdata$PAID_TIME))
table(is.na(Luminex.pdata$PAID_TIME))

table(Luminex.pdata$outcome,Luminex.pdata$timepoint)
####
##biomarker_pdata
PAMper.EC.marker<-read.csv(file = "rawdata/Data syndecan + thrombomodulin_PAMper.csv")
PAMper.Adiponectin.marker<-read.csv(file = "rawdata/data_Adiponectin_suPAR_Cell death_s100A10_PAMper.csv")
PAMper.VEGF<-read.csv(file = "rawdata/data_VEGF_R1_PAMper.csv")
PAMPer.biomarker<-inner_join(PAMper.EC.marker,PAMper.Adiponectin.marker,by="ï..Sample.ID") %>% inner_join(PAMper.VEGF,by="ï..Sample.ID")
PAMPer.biomarker<-PAMPer.biomarker[!duplicated(PAMPer.biomarker$ï..Sample.ID),]
timepoint<-PAMPer.biomarker$ï..Sample.ID %>% str_extract("(?<=\\-).+?(?=\\-)")
timepoint<- timepoint %>% gsub(pattern = "\\b0H\\b",replacement = "0HR") %>% gsub(pattern = "hr",replacement = "HR") %>% gsub(pattern = "Hr",replacement = "HR") %>% gsub(pattern = "27",replacement = "24") %>% gsub(pattern = "HR",replacement = "")
PAID<-PAMPer.biomarker$ï..Sample.ID %>% str_extract(".+?\\-") %>% str_extract(".+?\\-") %>% str_sub(5,-2)
PAMPer.biomarker<-PAMPer.biomarker[,-c(3,5,7,9,11,13,15,16,17)]

biomarker.pdata<-data.frame(rawname=PAMPer.biomarker$ï..Sample.ID,PAID=PAID, timepoint=timepoint)
biomarker.pdata$PAID_TIME<-paste(biomarker.pdata$PAID,biomarker.pdata$timepoint,sep = "_")
biomarker.pdata<-filter(biomarker.pdata,!is.na(timepoint))
biomarker.pdata$outcome<-Pamper.pdata$outcome[match(biomarker.pdata$PAID,Pamper.pdata$PAMPID)]
biomarker.pdata$outcome_time<-paste(biomarker.pdata$outcome,biomarker.pdata$timepoint,sep = "_")
dim(biomarker.pdata)
Filter_ID<-names(which(table(as.character(biomarker.pdata[which(biomarker.pdata$outcome %in% c("Non-resolving","Resolving")),"PAID"]))<2))

biomarker.pdata<-biomarker.pdata %>% filter(!outcome_time %in% c("Early-Nonsurvivors_24","Early-Nonsurvivors_72"))  %>% filter(!PAID %in% Filter_ID) 

table(duplicated(biomarker.pdata$PAID_TIME))
table(is.na(biomarker.pdata$PAID_TIME))

table(biomarker.pdata$outcome,biomarker.pdata$timepoint)

###Common_pdata
Common_PAID_TIME<-data_merge$PAID_TIME %>% intersect(Luminex.pdata$PAID_TIME) %>% intersect(biomarker.pdata$PAID_TIME)
Common_PAID<-Common_PAID_TIME %>% str_extract(".+?\\_") %>% str_sub(1,-2) %>% unique()

###Luminex exploration
Luminex.pdata<-Luminex.pdata[which(Luminex.pdata$PAID %in% Common_PAID ),]
rownames(Luminex.pdata)<-Luminex.pdata$rownames
PAMPer.luminex<-PAMPer.luminex[which(rownames(PAMPer.luminex) %in% Luminex.pdata$rownames),]
outcome_level<-c("Resolving_0","Non-resolving_0","Early-Nonsurvivors_0","Resolving_24","Non-resolving_24","Resolving_72","Non-resolving_72")

Luminex <- CreateSeuratObject(counts = t(PAMPer.luminex),meta.data = Luminex.pdata)
Luminex <- FindVariableFeatures(Luminex, selection.method = "vst", nfeatures = nrow(Luminex))
Luminex <- ScaleData(Luminex,do.scale = T,do.center = T)
###
###heatmap
Luminex$outcome_time<-factor(Luminex$outcome_time,levels = outcome_level)
Luminex.heatmap<-cbind(as.data.frame(scale(t(log2(as.data.frame(Luminex@assays$RNA@data)+1))))
                            ,outcome_timepoint=Luminex$outcome_time)
Luminex.heatmap2<-Luminex.heatmap %>% group_by(outcome_timepoint) %>% summarise_each(funs(mean)) 
timepoint<-Luminex.heatmap2$outcome_timepoint
Luminex.heatmap2<-Luminex.heatmap2[,-1]
Luminex.heatmap2<-t(Luminex.heatmap2)
colnames(Luminex.heatmap2)<-timepoint
col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))

ht<-Heatmap(Luminex.heatmap2,name = "Exp",col =col_fun,
            row_split = 2, show_row_names = T,cluster_columns = F,
            width = unit(10, "cm"),height = unit(10, "cm"))
ht


###biomarkers exploration
biomarker.pdata<-biomarker.pdata[which(biomarker.pdata$PAID %in% Common_PAID ),]
rownames(biomarker.pdata)<-biomarker.pdata$rawname
rownames(PAMPer.biomarker)<-PAMPer.biomarker$ï..Sample.ID
PAMPer.biomarker<-PAMPer.biomarker[,-1]
PAMPer.biomarker<-PAMPer.biomarker[which(rownames(PAMPer.biomarker) %in% biomarker.pdata$rawname),]
PAMPer.biomarker$Thrombomodulin_ng.ml<-as.numeric.factor(PAMPer.biomarker$Thrombomodulin_ng.ml)
PAMPer.biomarker$S100A10_ng.mL<-as.numeric.factor(PAMPer.biomarker$S100A10_ng.mL)
PAMPer.biomarker$suPAR_ng.mL<-as.numeric.factor(PAMPer.biomarker$suPAR_ng.mL)
PAMPer.biomarker$Cell.death_<-as.numeric.factor(PAMPer.biomarker$Cell.death_)
PAMPer.biomarker$Adiponectin_ng.mL<-as.numeric.factor(PAMPer.biomarker$Adiponectin_ng.mL)
PAMPer.biomarker$sFLT.1.sVEGFR1_pg.mL<-as.numeric.factor(PAMPer.biomarker$sFLT.1.sVEGFR1_pg.mL)
PAMPer.biomarker[is.na(PAMPer.biomarker)]<-0

outcome_level<-c("Resolving_0","Non-resolving_0","Early-Nonsurvivors_0","Resolving_24","Non-resolving_24")
biomarker <- CreateSeuratObject(counts = t(PAMPer.biomarker),meta.data = biomarker.pdata)
biomarker <- ScaleData(biomarker,do.scale = T,do.center = T,features = rownames(biomarker))

###
###heatmap
biomarker$outcome_time<-factor(biomarker$outcome_time,levels = outcome_level)
biomarker.heatmap<-cbind(as.data.frame(scale(t(log2(as.data.frame(biomarker@assays$RNA@data)+1))))
                       ,outcome_timepoint=biomarker$outcome_time)
biomarker.heatmap2<-biomarker.heatmap %>% group_by(outcome_timepoint) %>% summarise_each(funs(mean)) 
timepoint<-biomarker.heatmap2$outcome_timepoint
biomarker.heatmap2<-biomarker.heatmap2[,-1]
biomarker.heatmap2<-t(biomarker.heatmap2)
colnames(biomarker.heatmap2)<-timepoint
col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))

ht<-Heatmap(biomarker.heatmap2,name = "Exp",col =col_fun,
            row_split = 2, show_row_names = T,cluster_columns = F,
            width = unit(10, "cm"),height = unit(10, "cm"))
ht
####Visulization for individual cytokines
cytokine.data<-FetchData(Luminex,vars = c(rownames(Luminex),"outcome","timepoint","outcome_time"))
Luminex$outcome_time
####
cytokine.lineplot<-list()
cytokine.name<-list()
for(i in c(1:21)){  
  cytokine.name<-colnames(cytokine.data)[i]
  p.tmp<-ggline(cytokine.data, x = "timepoint", y = cytokine.name, 
              add = c("median_mad"),size = 1,point.size = 0.5,
              color = "outcome")+theme(aspect.ratio=1,legend.position = "bottom")
  cytokine.lineplot[[cytokine.name]]<-p.tmp
  }


ggarrange(plotlist = cytokine.lineplot[Subset3], 
          ncol = 3, nrow = 3,common.legend = T,align = "v")
####
cytokine.boxplot<-list()
cytokine.name<-list()
for(i in c(1:21)){  
  cytokine.name<-colnames(cytokine.data)[i]
  my_comparisons <- list( c("Resolving_0", "Non-resolving_0"), c("Resolving_0", "Early-Nonsurvivors_0"), c("Resolving_24", "Non-resolving_24"),c("Resolving_72","Non-resolving_72"))
  p.tmp<-ggboxplot(cytokine.data, x = "outcome_time", y = cytokine.name,
            color = "outcome_time", palette = "jco")+theme(aspect.ratio=1,legend.position = "bottom",axis.text.x =  element_text(angle = 90),axis.title.x=element_blank())+  stat_compare_means(comparisons = my_comparisons)+ scale_y_continuous(trans='log10')
  cytokine.boxplot[[cytokine.name]]<-p.tmp
}


ggarrange(plotlist = cytokine.boxplot[Subset3], 
          ncol = 4, nrow = 2,common.legend = T,align = "v")
?ggboxplot
 