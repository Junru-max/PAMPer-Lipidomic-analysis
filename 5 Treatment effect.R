library(imputeMissings)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(dendsort)
library(tidyverse)
library(ggpubr)
library(rstatix)
####Fig3 DEFG

####Data manipulation
table( pdata.trauma$PAMPID %in% sum.data$PAID )
match(pdata.trauma$PAMPID,sum.data$PAID)
sum.data.0h<-sum.data[which(sum.data$TIME=="0"),]

pdata.trauma$early_death<-0
pdata.trauma$early_death[which(pdata.trauma$outcome=="Early-Nonsurvivors")]<-1
pdata.trauma$early_death<-factor(pdata.trauma$early_death,levels = c(0,1))
pdata.trauma$Total_SFA<-sum.data.0h$Total_SFA[match(pdata.trauma$PAMPID,sum.data.0h$PAID)]
pdata.trauma$Total_USFA<-sum.data.0h$Total_USFA[match(pdata.trauma$PAMPID,sum.data.0h$PAID)]
pdata.trauma$Total_FA<-pdata.trauma$Total_SFA+pdata.trauma$Total_USFA
pdata.trauma$PE_USFA<-sum.data.0h$PE_USFA[match(pdata.trauma$PAMPID,sum.data.0h$PAID)]
pdata.trauma$RR_IMP<-pdata.trauma$vitals_resp_rate
pdata.trauma$RR_IMP[is.na(pdata.trauma$RR_IMP)]<-median(pdata.trauma$RR_IMP[!is.na(pdata.trauma$RR_IMP)])
pdata.trauma$age_index<-cut(pdata.trauma$age, breaks = c(0,55,max(pdata.trauma$age)),labels = 0:1)
pdata.trauma$age_index<-as.numeric.factor(pdata.trauma$age_index)

pdata.trauma$GCS_RTS<-cut(pdata.trauma$initial_GCS, breaks = c(0,3,5,8,12,15),labels = 0:4)
pdata.trauma$GCS_RTS<-as.numeric.factor(pdata.trauma$GCS_RTS)
pdata.trauma$SBP_RTS<-cut(pdata.trauma$vitals_sbp, breaks = c(-1,0,49,75,89,max(pdata.trauma$vitals_sbp)),labels = 0:4)
pdata.trauma$SBP_RTS<-as.numeric.factor(pdata.trauma$SBP_RTS)
pdata.trauma$RR_RTS<-cut(pdata.trauma$RR_IMP, breaks = c(-1,0,5,9,29,max(pdata.trauma$RR_IMP)),labels = c(0,1,2,4,3))
pdata.trauma$RR_RTS<-as.numeric.factor(pdata.trauma$RR_RTS)
pdata.trauma$RTS<- 0.9368*pdata.trauma$GCS_RTS + 0.7326 *pdata.trauma$GCS_RTS + 0.2908 * pdata.trauma$RR_RTS
hist(pdata.trauma$RTS,breaks = 20)

pdata.trauma$b_TRISS_blunt<- 0.8085*pdata.trauma$RTS-0.0835*pdata.trauma$iss- 1.743*pdata.trauma$age_index -0.4499
pdata.trauma$b_TRISS_pene<- 0.9934*pdata.trauma$RTS- 0.0651*pdata.trauma$iss- 1.1360 *pdata.trauma$age_index -2.5355
pdata.trauma$b_TRISS<-pdata.trauma$b_TRISS_pene
pdata.trauma$b_TRISS[which(pdata.trauma$injury_type=="1:Blunt")]<-pdata.trauma$b_TRISS_blunt[which(pdata.trauma$injury_type=="1:Blunt")]
pdata.trauma$mortality_TRISS<- 1- (1/(1+ exp(-pdata.trauma$b_TRISS)))
pdata.trauma$injury_type
hist(pdata.trauma$mortality_TRISS,breaks = 100)
pdata.trauma
table(pdata.trauma$age_index,pdata.trauma$age)

pdata.trauma$early_death
####Exploration 
library(ggplot2)
pdata.trauma$Total_FA_CA<-cut(pdata.trauma$Total_FA, breaks = c(0,median(pdata.trauma$Total_FA),max(pdata.trauma$Total_FA)),labels = c("Low","High"))
pdata.trauma$Total_FA_CA<-relevel(pdata.trauma$Total_FA_CA,ref = "High")
pdata.trauma$TREATMENT<-factor(pdata.trauma$TREATMENT,levels = c("standard care","prehospital plasma"))
cor(pdata.trauma$mortality_TRISS,pdata.trauma$Total_FA,method = "spearman")
ggplot(pdata.trauma, aes(x=mortality_TRISS, y=Total_FA, color=TREATMENT, shape=TREATMENT)) +
  geom_point() + 
  geom_smooth(method=loess, se=FALSE, fullrange=TRUE)+ ylim(0,2500)+ geom_text(aes(label=ifelse(early_death==1,"D",'')),hjust=0,vjust=0)+geom_hline(yintercept=median(pdata.trauma$Total_FA), linetype="dashed", 
                                                                                                                                                    color = "black", size=1)+geom_vline(xintercept=0.5, linetype="dashed",                                                                                                                                                                                       color = "black", size=1)+theme(aspect.ratio=1,legend.position = "bottom")
####
table(pdata.trauma$early_death[which(pdata.trauma$mortality_TRISS<0.5 & pdata.trauma$Total_FA>median(pdata.trauma$Total_FA) )])



pdata.trauma$traumatic_brain_injury<-PAMPer.imp$`TRAUMATIC BRAIN INJURY`[match(pdata.trauma$PAMPID,PAMPer.imp$`PAMPER ID NUMBER`)] %>% droplevels()
table(pdata.trauma$traumatic_brain_injury,pdata.trauma$severe_head)



#####Logistic regression
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(geepack)
model1<-glm(early_death~gender+traumatic_brain_injury+mortality_TRISS+Total_FA_CA+TREATMENT,data=pdata.trauma,family="binomial" )
summary(model1)
plot_model(model1, transform = NULL,value.offset = .3,show.values = TRUE)+theme(aspect.ratio = 1.5)
####GEE for Logistic regression
pdata.trauma2<-select(pdata.trauma, early_death, gender, traumatic_brain_injury, mortality_TRISS, Total_FA_CA,
                       SiteID,base, TREATMENT)
pdata.trauma2 <- na.omit(pdata.trauma2)
pdata.trauma2$early_death<-as.numeric.factor(pdata.trauma2$early_death)
str(pdata.trauma2)
table(pdata.trauma2$early_death)
model2<-geeglm(early_death~gender+traumatic_brain_injury+mortality_TRISS+Total_FA_CA+TREATMENT,data=pdata.trauma2,family=binomial,id = SiteID,corstr = "exchangeable" )

summary(model2)
plot_model(model2, transform = NULL,value.offset = .3,show.values = TRUE)+theme(aspect.ratio = 1.5)
model3<-geeglm(early_death~gender+traumatic_brain_injury+mortality_TRISS+Total_FA_CA+TREATMENT,data=pdata.trauma2,family=binomial,id = base,corstr = "exchangeable" )
plot_model(model3, transform = NULL,value.offset = .3,show.values = TRUE)+theme(aspect.ratio = 1.5)



####Correlation matrix
table(pdata.trauma$initial_GCS,pdata.trauma$initial_GCS_8)
####For casual inference of early-outcome
library(Seurat)
pdata.trauma$Total_log2FA<-scale(log2(pdata.trauma$Total_FA+1))
pdata.trauma$INR[is.na(pdata.trauma$INR)]<-median(pdata.trauma$INR[!is.na(pdata.trauma$INR)])
pdata.trauma$PAID<-as.character(pdata.trauma$PAMPID)
pdata.trauma$Total_log2FA
#Varabile<-c("TREATMENT","age","iss","initial_GCS_8","early_death","PAID","Total_FA_CA","mortality_TRISS","traumatic_brain_injury","ed_coagulopathy","PH_crystalloid","transfusion_24h","crystalloid_24h","PH_intubation","INR")
Varabile2<-c("TREATMENT","age","iss","initial_GCS","early_death","PAID","Total_log2FA","mortality_TRISS","traumatic_brain_injury","ed_coagulopathy","PH_crystalloid","transfusion_24h","crystalloid_24h","PH_intubation","INR")

raw.beya.matrix<-pdata.trauma[,Varabile2]
Idents(biomarker)<-biomarker$timepoint
biomarker.0h.heatmap<-FetchData(subset(biomarker,idents="0"),vars = c(rownames(biomarker)))
biomarker.0h.heatmap[is.na(biomarker.0h.heatmap)]<- 0
biomarker.0h.heatmap<-as.data.frame(scale(log2(biomarker.0h.heatmap+1)))
biomarker.0h.heatmap<-cbind(biomarker.0h.heatmap,PAID=as.character(subset(biomarker,idents="0")$PAID))
Idents(Luminex)<-Luminex$timepoint
Luminex.heatmap<-cbind(as.data.frame(scale(t(log2(as.data.frame(subset(Luminex,idents="0")@assays$RNA@data)+1))))
                       ,PAID=as.character(subset(Luminex,idents="0")$PAID))

raw.beya.matrix<-biomarker.0h.heatmap %>% left_join(raw.beya.matrix,by="PAID") %>% left_join(Luminex.heatmap,by="PAID")
raw.beya.matrix<-raw.beya.matrix[,-which(colnames(raw.beya.matrix)=="PAID")]
raw.beya.matrix$TREATMENT<-factor(as.character(raw.beya.matrix$TREATMENT),levels = c("prehospital plasma","standard care"),labels = c(1,2))
raw.beya.matrix$TREATMENT<-fastDummies::dummy_cols(raw.beya.matrix$TREATMENT)[,3]
raw.beya.matrix$early_death<-fastDummies::dummy_cols(raw.beya.matrix$early_death)[,3]
raw.beya.matrix$traumatic_brain_injury<-fastDummies::dummy_cols(raw.beya.matrix$traumatic_brain_injury)[,3]
raw.beya.matrix$Total_FA_CA<-fastDummies::dummy_cols(raw.beya.matrix$Total_FA_CA)[,3]
raw.beya.matrix$PH_intubation<-fastDummies::dummy_cols(raw.beya.matrix$PH_intubation)[,3]
raw.beya.matrix$ed_coagulopathy<-fastDummies::dummy_cols(raw.beya.matrix$ed_coagulopathy)[,3]
raw.beya.matrix$initial_GCS_8<-fastDummies::dummy_cols(raw.beya.matrix$initial_GCS_8)[,3]
saveRDS(raw.beya.matrix,file = "process_PAMPer.Rdata")
###correlation matrix
cor_matrix<-cor(raw.beya.matrix,method = "spearman")
cor_matrix[is.na(cor_matrix)]<-0
cor_matrix<-as.data.frame(cor_matrix)

Subset1<-c("IL.6","IL.8","IL.10","MCP.1","MIG","IP.10","TNFa")
Subset2<-c("IL.4","IL.5","IL.1b", "IL.2","IL.7", "IL.17A","GM.CSF")
Subset3<-c("IL.17E","IL.22","IL.9","IL.33","IL.21","IL.23.ng.ml.","IL.27")
EC<-c("Cell.death-","S100A10-ng.mL","sFLT.1.sVEGFR1-pg.mL","Syndecan.1-ng.ml","Thrombomodulin-ng.ml","suPAR-ng.mL")
Adipokine<-c("Adiponectin-ng.mL")
ARMS<-"TREATMENT"
Injury_Severity<-c("iss","traumatic_brain_injury","initial_GCS_8","mortality_TRISS")
Coagulation<-c("ed_coagulopathy","INR")
Other_Intervention<-c("PH_intubation","transfusion_24h","crystalloid_24h")
Age<-"age"
Lipid<-"Total_FA_CA"
outcome<-"early_death"
group<-rownames(cor_matrix)
group[which(group %in% Subset1)]<-"Subset1 Cytokines"
group[which(group %in% Subset2)]<-"Subset2 Cytokines"
group[which(group %in% Subset3)]<-"Subset3 Cytokines"
group[which(group %in% EC)]<-"EC injury biomarkers"
group[which(group %in% Adipokine)]<-"Adipokine"
group[which(group %in% Subset3)]<-"Subset3 Cytokines"
group[which(group %in% ARMS)]<-"Treatment Arms"
group[which(group %in% Injury_Severity)]<-"Injury Severity"
group[which(group %in% Coagulation)]<-"Coagulation"
group[which(group %in% Other_Intervention)]<-"Other Intervention"
group[which(group %in% Age)]<-"Age"
group[which(group %in% Lipid)]<-"Lipid"
group[which(group %in% outcome)]<-"Outcome"
group<-factor(group)
myPalette2<-hue_pal()(length(unique(group)))
names(myPalette2)<-levels(group)
ha = rowAnnotation(Class= group,
                   col = list(Class=myPalette2))

col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))

ht<-Heatmap(cor_matrix,name = "Exp",col =col_fun, show_row_names = T,cluster_columns = T,
            width = unit(15, "cm"),height = unit(15, "cm"),right_annotation = ha,show_row_dend = F,show_column_dend = T)
ht
ht<-Heatmap(cor_matrix,name = "Exp",col =col_fun,cluster_rows = T,
            show_row_names = T,cluster_columns = T,width = unit(20, "cm"),height = unit(20, "cm"),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.1f", cor_matrix[i, j]), x, y, gp = gpar(fontsize = 5))
            })
ht
table(raw.beya.matrix$Total_FA_CA,raw.beya.matrix$TREATMENT
)

sort(cor_matrix$Total_FA_CA)
#####Data preparation for causual inference

normality<-shapiro_test(raw.beya.matrix,vars = colnames(raw.beya.matrix))
hist(normality$p,breaks = 20)
ggqqplot(raw.beya.matrix,colnames(raw.beya.matrix)[c(1:7,13,24:40)],combine = T)

colnames(raw.beya.matrix)<-c("Syndecan-1","Thrombomodulin","S100A10","suPAR","HC-DNA","Adiponectin","sVEGFR1","FFP","Age","ISS","GCS","Early-Death","Total_Lipid_Conc.",
                             "Mortality_TRISS","TBI","Coagulopathy","PH_crystalloid","Transfusion_24h","Crystalloid_24h","PH_intubation","INR","GM-CSF","IL-10","IL-17A","IL-1B",
                             "IL-2","IL-4","IL-5","IL-6","IL-7","IL-8","IP-10","MCP-1","TNFa","IL-22","IL-9","IL-33","IL-21","IL-23","IL-17E","IL-27","MIG")

write.table(raw.beya.matrix,file = "PAMPer.txt",sep = "\t",row.names = F,col.names = T,quote = F)
#####Data preparation for causual inference2
Subset1<-c("IL-6","IL-8","IL-10","MCP-1","MIG","IP-10","TNFa")
Subset2<-c("IL-4","IL-5","IL-1B", "IL-2","IL-7", "IL-17A","GM-CSF")
Subset3<-c("IL-17E","IL-22","IL-9","IL-33","IL-21","IL-23","IL-27")
raw.beya.matrix$CS1<-apply(raw.beya.matrix[,Subset1],1,mean)
raw.beya.matrix$CS2<-apply(raw.beya.matrix[,Subset2],1,mean)
raw.beya.matrix$CS3<-apply(raw.beya.matrix[,Subset3],1,mean)
raw.beya.matrix<-raw.beya.matrix[,-which(colnames(raw.beya.matrix) %in% c(Subset1,Subset2,Subset3))]

write.table(raw.beya.matrix,file = "PAMPer.txt",sep = "\t",row.names = F,col.names = T,quote = F)

####therotical total lipid conc. with actually count
Conc.median<-median(PAMPer.imp$Total.conc[which(PAMPer.imp$outcome=="Baseline")])

#pdata.trauma$weight[which(pdata.trauma$weight_units=="1:Pounds")]<-0.45359237*pdata.trauma$weight[which(pdata.trauma$weight_units=="1:Pounds")]
#pdata.trauma$weight[is.na(pdata.trauma$weight)]<-median(pdata.trauma$weight[!is.na(pdata.trauma$weight)])
pdata.trauma$Blood_volume<-0.075*pdata.trauma$weight
pdata.trauma$Blood_volume[which(pdata.trauma$gender=="2:Female")]<-pdata.trauma$Blood_volume[which(pdata.trauma$gender=="2:Female")]-0.5
pdata.trauma$Blood_volume[which(pdata.trauma$Blood_volume>10)]<-10
hist(pdata.trauma$Blood_volume,breaks = 20)
pdata.trauma$Total_FA_T<-pdata.trauma$Total_FA
pdata.trauma$Total_FA_T[which(pdata.trauma$TREATMENT=="prehospital plasma")]<-pdata.trauma$Total_FA_T[which(pdata.trauma$TREATMENT=="prehospital plasma")]-0.5*Conc.median/pdata.trauma$Blood_volume[which(pdata.trauma$TREATMENT=="prehospital plasma")]
wilcox.test(Total_FA_T~TREATMENT,data = pdata.trauma[which(pdata.trauma$outcome!="Early-Nonsurvivors"),])
wilcox.test(Total_FA~TREATMENT,data = pdata.trauma[which(pdata.trauma$outcome!="Early-Nonsurvivors"),])
###
Stat.casual$Number.of.directed.Edges.
Stat.casual<-read.csv(file = "rawdata/Statisictis for causal inference.csv")
Stat.casual<-Stat.casual %>% pivot_longer(
  cols = Number.of.nodes:Number.of.directed.Edges.,
  names_to = "Types",
  values_to = "count"
)
ggline((Stat.casual %>% filter(r2==0.2 & r3 ==0.2 & ï..Alpha==0.05 )), x = "r1", y = "count",size = 1,point.size = 3, color = "Types")+theme(aspect.ratio=1,legend.position = "bottom")
ggline((Stat.casual %>% filter(r1==0.2 & r3 ==0.2 & ï..Alpha==0.05 )), x = "r2", y = "count",size = 1,point.size = 3, color = "Types")+theme(aspect.ratio=1,legend.position = "bottom")
ggline((Stat.casual %>% filter(r1==0.2 & r2 ==0.2 & ï..Alpha==0.05 )), x = "r3", y = "count",size = 1,point.size = 3, color = "Types")+theme(aspect.ratio=1,legend.position = "bottom")

