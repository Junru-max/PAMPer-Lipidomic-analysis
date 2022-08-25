library(Seurat)
library(dplyr)
library(tidyr)
library(foreign)
library(stringr)
library(dplyr)
###MM
MM.pdata<-as.data.frame(t(read.csv(file = "Rawdata/pdata patients.csv",row.names = 1)))
MM.pdata$TIME_POINT<-as.character(MM.pdata$TIME_POINT)
MM.pdata$TIME_POINT[which(MM.pdata$TIME_POINT=="*D2")]<-"D2-D5"
MM.pdata$outcome<-"Non-resolving"
MM.pdata$outcome[which(MM.pdata$GROUP_ID %in% c("UI_MM","Inf_MM"))]<-"Resolving"
table(MM.pdata$GROUP_ID)
MM.pdata$GROUP.TIME<-paste(MM.pdata$outcome,MM.pdata$TIME_POINT,sep = "_")
MM.pdata$GROUP.TIME<-factor(MM.pdata$GROUP.TIME,levels = c("Resolving_<6hr","Non-resolving_<6hr","Resolving_24hr",
                                                       "Non-resolving_24hr","Resolving_D2-D5","Non-resolving_D2-D5"))
table(MM.pdata$TIME,MM.pdata$TIME_POINT)

MM.lipid.pdata<-read.csv(file = "Rawdata/MM.lipid.pdata.csv",row.names = 1)
MM.lipid.pdata$PAMPer_Lipid<-MM.lipid.pdata$BIOCHEMICAL
MM.lipid.pdata<-MM.lipid.pdata[!duplicated(MM.lipid.pdata$PAMPer_Lipid),]
MM.lipid.pdata<-MM.lipid.pdata[MM.lipid.pdata$PAMPer_Lipid %in% rownames(PAMPer.imp),]

MM.lipid.pdata$SUB_PATHWAY<-as.character(MM.lipid.pdata$SUB_PATHWAY)
MM.lipid.pdata$SUB_PATHWAY<-factor(MM.lipid.pdata$SUB_PATHWAY,levels = c("Monoacylglycerol","Diacylglycerol","Phosphatidylethanolamine (PE)","Phosphatidylcholine (PC)",
                                                                     "Phosphatidylinositol (PI)","Ceramides","Hexosylceramides (HCER)","Lactosylceramides (LCER)","Sphingomyelins"))
###
MM.impdata<-read.csv(file = "Rawdata/VolNormImpData.csv",header = T,row.names = 1)
MM.impdata<-MM.impdata[rownames(MM.lipid.pdata),]
rownames(MM.impdata)<-MM.lipid.pdata$PAMPer_Lipid[match(rownames(MM.impdata),rownames(MM.lipid.pdata))]
MM.imp<-CreateSeuratObject(counts = MM.impdata)
MM.imp<-AddMetaData(MM.imp,MM.pdata)
MM.imp <- ScaleData(MM.imp,do.scale = T,do.center = T,features = rownames(MM.imp))

###Covid-1
##pdata
covid.pdata<-read.csv(file = "rawdata/covid_pdata.csv")
covid.pdata$TIMEPOINT<-"N"
covid.pdata$TIMEPOINT[which(covid.pdata$Days.after.progression <3 & covid.pdata$Group.d==2)]<-"<48h"
covid.pdata$TIMEPOINT[which(covid.pdata$Days.after.progression >2 & covid.pdata$Days.after.progression<6 & covid.pdata$Group.d==2)]<-"D2-D5"
covid.pdata$TIMEPOINT[which(covid.pdata$Days.after.progression >5 & covid.pdata$Group.d==2)]<-"D6_D14"
covid.pdata$TIMEPOINT[which(covid.pdata$Days.after.progression <0 & covid.pdata$Group.d==3)]<-"Before"
covid.pdata$TIMEPOINT[which(covid.pdata$Days.after.progression > -1 & covid.pdata$Days.after.progression <3 & covid.pdata$Group.d==3)]<-"<48h"
covid.pdata$TIMEPOINT[which(covid.pdata$Days.after.progression >  2 & covid.pdata$Days.after.progression <6 & covid.pdata$Group.d==3)]<-"D2-D5"
covid.pdata$TIMEPOINT[which(covid.pdata$Days.after.progression >  5 & covid.pdata$Group.d==3)]<-"D6_D14"
covid.pdata$TIMEPOINT[which(covid.pdata$Group.d==0)]<-"Baseline"
covid.pdata$TIMEPOINT[which(covid.pdata$Group.d==1)]<-"Non-Covid"
table(covid.pdata$TIMEPOINT)
covid.pdata$outcome<-covid.pdata$Group.d
covid.pdata$GROUP.TIME<-paste(covid.pdata$Group.d,covid.pdata$TIMEPOINT,sep = "_")
covid.pdata$GROUP.TIME2<-paste(covid.pdata$Group.d,covid.pdata$Days.after.progression,sep = "_")

covid.pdata$GROUP.TIME<-factor(covid.pdata$GROUP.TIME,levels = c("0_Baseline","1_Non-Covid","2_<48h","2_D2-D5","2_D6_D14","3_Before","3_<48h","3_D6_D14"),
                             labels = c("Baseline","Non-covid","Non-severe_<48h","Non-severe_D2-D5","Non-severe_D6-D14","Severe_before_pro","Severe_after_<48h","Severe_after_D6-D14"))
covid.pdata$outcome<-factor(covid.pdata$outcome,levels = c(0,1,2,3),labels = c("Baseline","Non-covid","Mild-Mod Covid","Severe Covid"))
covid.pdata$TIME<-factor(covid.pdata$TIME,levels = c("Baseline","Non-Covid","Before","<48h","D2-D5","D6_D14"))

##concentration matrix
covid.data<-readRDS(file = "rdata/Covid19.log2.data.Rdata")
covid.data[1:5,1:5]
table(rownames(Com.lipid.list) %in% rownames(covid.data))
covid.pdata<-covid.pdata[which(covid.pdata$Metabolomics.ID.e %in% colnames(covid.data)),]

covid.lipid.pdata<-MM.lipid.pdata[which(rownames(MM.lipid.pdata) %in% rownames(covid.data)),]

covid.data<-covid.data[rownames(covid.lipid.pdata),as.character(covid.pdata$Metabolomics.ID.e)]
rownames(covid.pdata)<-covid.pdata$Metabolomics.ID.e
rownames(covid.data)<-covid.lipid.pdata$PAMPer_Lipid

covid.imp<-CreateSeuratObject(counts = covid.data) %>% AddMetaData(covid.pdata)
Idents(covid.imp)<-covid.imp$GROUP.TIME2
covid.imp<-subset(covid.imp,idents="3_9",invert=T)
covid.imp<-ScaleData(covid.imp,features = rownames(covid.imp))


###Covid-2
##pdata
covid.2.plasma<-read.csv(file = "rawdata/covid_2_plasma.csv",row.names = 1)
covid.2.plasma<-log2(covid.2.plasma+1)
covid.2.plasma<-covid.2.plasma[-grep("p",rownames(covid.2.plasma)),]
covid.2.pdata<-read.csv(file = "rawdata/covid_2_pdata.csv")

covid.2.plasma.TAG<-rownames(covid.2.plasma)[grep("TAG",rownames(covid.2.plasma))]
covid.2.plasma.TAG<-covid.2.plasma.TAG[-221]
covid.2.plasma.TAG2<-covid.2.plasma.TAG %>% str_sub(0,-7) 
covid.2.plasma.TAG3<-covid.2.plasma.TAG %>% str_sub(9,-2) 
covid.2.plasma.TAG4<-paste(covid.2.plasma.TAG2,covid.2.plasma.TAG3,sep = "-FA")

covid.2.plasma.DAG<-rownames(covid.2.plasma)[grep("DAG",rownames(covid.2.plasma))]
covid.2.plasma.DAG<-covid.2.plasma.DAG[-29]
covid.2.plasma.DAG2<-covid.2.plasma.DAG %>% str_sub(0,-16) 
covid.2.plasma.DAG3<-covid.2.plasma.DAG %>% str_sub(8,-1)
covid.2.plasma.DAG4<-paste(covid.2.plasma.DAG2,covid.2.plasma.DAG3,sep = "")

covid.2.plasma.PE<-rownames(covid.2.plasma)[grep("^PE",rownames(covid.2.plasma))]
covid.2.plasma.PE<-covid.2.plasma.PE[1:23]
covid.2.plasma.PE2<-covid.2.plasma.PE %>% str_sub(0,-16) 
covid.2.plasma.PE3<-covid.2.plasma.PE %>% str_sub(7,-1)
covid.2.plasma.PE4<-paste(covid.2.plasma.PE2,covid.2.plasma.PE3,sep = "")

covid.2.plasma.PC<-rownames(covid.2.plasma)[grep("^PC",rownames(covid.2.plasma))]
covid.2.plasma.PC<-covid.2.plasma.PC[1:43]
covid.2.plasma.PC2<-covid.2.plasma.PC %>% str_sub(0,-16) 
covid.2.plasma.PC3<-covid.2.plasma.PC %>% str_sub(7,-1)
covid.2.plasma.PC4<-paste(covid.2.plasma.PC2,covid.2.plasma.PC3,sep = "")

covid.2.plasma.PI<-rownames(covid.2.plasma)[grep("^PI",rownames(covid.2.plasma))]
covid.2.plasma.PI<-covid.2.plasma.PI[1:14]
covid.2.plasma.PI2<-covid.2.plasma.PI %>% str_sub(0,-17) 
covid.2.plasma.PI3<-covid.2.plasma.PI %>% str_sub(8,-1)
covid.2.plasma.PI4<-paste(covid.2.plasma.PI2,covid.2.plasma.PI3,sep = "")

covid.2.plasma.SM<-rownames(covid.2.plasma)[grep("^SM",rownames(covid.2.plasma))]
covid.2.plasma.SM<-covid.2.plasma.SM[grep("d18:1",covid.2.plasma.SM)]
covid.2.plasma.SM4<-covid.2.plasma.SM %>% gsub(pattern = " d18:1/",replacement = "(")
covid.2.plasma.SM4<-paste(covid.2.plasma.SM4,rep(")",length(covid.2.plasma.SM4)),sep = "")


covid2.lipid.pdata<-data.frame(covid2id=c(covid.2.plasma.TAG,covid.2.plasma.DAG,covid.2.plasma.PE,covid.2.plasma.PC,covid.2.plasma.PI,covid.2.plasma.SM),
                                   PAMPer.ID=c(covid.2.plasma.TAG4,covid.2.plasma.DAG4,covid.2.plasma.PE4,covid.2.plasma.PC4,covid.2.plasma.PI4,covid.2.plasma.SM4)
)
covid2.lipid.pdata<-covid2.lipid.pdata[!duplicated(covid2.lipid.pdata$PAMPer.ID),]
covid2.lipid.pdata<-covid2.lipid.pdata[which(covid2.lipid.pdata$PAMPer.ID %in% covid.lipid.pdata$PAMPer_Lipid),]

covid.2.plasma<-covid.2.plasma[as.character(covid2.lipid.pdata$covid2id),]
rownames(covid.2.plasma)<-covid2.lipid.pdata$PAMPer.ID
covid2.imp<-CreateSeuratObject(counts = covid.2.plasma)
covid2.imp$GROUP.TIME<-covid.2.pdata$GroupID[match(colnames(covid.2.plasma),covid.2.pdata$ï..SampleID)]
table(covid2.imp$GROUP.TIME)
covid2.imp <- ScaleData(covid2.imp,do.scale = T,do.center = T,features = rownames(covid2.imp))
####
Com.lipid.list<-MM.lipid.pdata[which(MM.lipid.pdata$BIOCHEMICAL %in% rownames(covid2.imp)),]
Com.lipid.list$Covid.2.lipid<-covid2.lipid.pdata$covid2id[match(Com.lipid.list$PAMPer_Lipid,covid2.lipid.pdata$PAMPer.ID)]
write.csv(Com.lipid.list,file = "figure/20210421 Figures for reviewers/Com.lipid.list.csv")



###Covid-3
##pdata
