library(Seurat)
library(dplyr)
library(ggpubr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(dendsort)
library(stringr)
library(tidyverse)
library(rstatix)
library(emmeans)
###Fig2 D.E.F
######Lipid_class
data_class<-read.csv(file = "rawdata/data_lipidclass.csv",header = T,row.names = 1)
data_class[1:5,1:5]
data_class<-data_class[,pdata$ID]
data_class<-CreateSeuratObject(counts = data_class,meta.data = pdata)
data_class <- FindVariableFeatures(data_class, selection.method = "vst", nfeatures = nrow(data_class))
data_class <- ScaleData(data_class,do.scale = T,do.center = T)
###heatmap
data_class$outcome_time
PAMPer.heatmap<-FetchData(data_class,slot = "scale.data",vars = c(rownames(data_class),"outcome_time"))
PAMPer.heatmap2<-PAMPer.heatmap %>% group_by(outcome_time) %>% summarise_each(funs(mean)) 
timepoint<-PAMPer.heatmap2$outcome_time
PAMPer.heatmap2<-PAMPer.heatmap2[,-1]
PAMPer.heatmap2<-t(PAMPer.heatmap2)
colnames(PAMPer.heatmap2)<-timepoint

col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))

ht<-Heatmap(PAMPer.heatmap2,name = "Exp",col =col_fun,
             show_row_names = T,cluster_columns = F,cluster_rows = T,width = unit(10, "cm"),height = unit(10, "cm"))
ht
####Lipid_species
data_species<-read.csv(file = "rawdata/data_fattyacid.csv",header = T,row.names = 1)
data_species<-data_species[,pdata$ID]
data_species<-CreateSeuratObject(counts = data_species,meta.data = pdata)
data_species <- ScaleData(data_species,do.scale = T,do.center = T)

###heatmap

total_fd<-rownames(data_species)[grep("Total",rownames(data_species))]
PAMPer.heatmap<-FetchData(data_species,slot = "scale.data",vars = c(total_fd,"outcome_time"))
PAMPer.heatmap2<-PAMPer.heatmap %>% group_by(outcome_time) %>% summarise_each(funs(mean)) 
timepoint<-PAMPer.heatmap2$outcome_time
PAMPer.heatmap2<-PAMPer.heatmap2[,-1]
PAMPer.heatmap2<-t(PAMPer.heatmap2)
colnames(PAMPer.heatmap2)<-timepoint

col_fun = colorRamp2(c(-1, 0,1), c("slateblue2", "white", "firebrick2"))

ht<-Heatmap(PAMPer.heatmap2,name = "Exp",col =col_fun,
            show_row_names = T,cluster_columns = F)
ht
###

####summarized data
sum.data<-read.csv(file = "rawdata/summarized_data.csv",row.names = 1)
sum.data<-sum.data[,pdata$ID]
sum.data[is.na(sum.data)]<-0
species=rownames(sum.data) %>% str_extract(".+?\\[")  %>% 
  str_extract(".+?\\[") %>% str_sub(1,-2)
Fattyacid=rownames(sum.data) %>% str_extract("(?<=\\[).+?(?=\\])")
sum.pdata<-data.frame(species=species,Fattyacid=Fattyacid )
sum.pdata$Fa.group<-"SFA"
sum.pdata$Fa.group[which(sum.pdata$Fattyacid %>% str_sub(6,-1) !="0")]<-"USFA"
rownames(sum.pdata)<-rownames(sum.data)
sum.pdata$species.group<-paste(sum.pdata$species,sum.pdata$Fa.group,sep = "_")
sum.data$species.group<-sum.pdata$species.group

sum.data<-sum.data %>% group_by(species.group) %>% summarise_each(funs(sum))
species.group<-sum.data$species.group
sum.data<-sum.data[,which(colnames(sum.data)!="species.group")]
sum.data<-as.data.frame(t(sum.data))
colnames(sum.data)<-species.group
sum.data$outcome_time<-PAMPer.imp$outcome_time
sum.data$TREATMENT<-PAMPer.imp$TREATMENT
sum.data$outcome<-PAMPer.imp$outcome
sum.data$TIME<-PAMPer.imp$`TIME POINT`
sum.data$PAID<-PAMPer.imp$`PAMPER ID NUMBER`
ggline(sum.data, x = "TIME", y = "PE_SFA", 
       add = c("median_mad","jitter"),size = 1,point.size = 3,
       color = "outcome")+theme(aspect.ratio=1,legend.position = "bottom")
sum.data$PE_USFA
#######statistical analysis For outcome_time 

sum.data.ex.nons<-sum.data[which(sum.data$outcome %in% c("Resolving","Non-resolving")),] %>% droplevels()
sum.data.ex.nons %>%
  group_by(outcome,TIME) %>%
  shapiro_test(DAG_USFA)
ggqqplot(sum.data.ex.nons, "DAG_USFA", ggtheme = theme_bw()) +
  facet_grid(TIME ~ outcome, labeller = "label_both")
sum.data.ex.nons<-tibble(sum.data.ex.nons)
res.aov <-sum.data.ex.nons %>% anova_test(
  DAG_USFA ~ outcome*TIME )
#####2-way anova
get_anova_table(res.aov)

model <- lm(DAG_USFA ~ TIME * outcome, data = sum.data.ex.nons)
sum.data.ex.nons %>%
  group_by(TIME) %>%
  anova_test(DAG_USFA ~ outcome, error = model)
##
# pairwise comparisons
pwc <- sum.data.ex.nons %>% 
  group_by(TIME) %>%
  emmeans_test(DAG_USFA ~ outcome, p.adjust.method = "BH") 
pwc

####kruskal.test
sum.data.ex.nons<-sum.data[which(sum.data$outcome_time %in% c("Resolving_0","Non-resolving_0","Early-Nonsurvivors_0")),] %>% droplevels()

kruskal.test(DAG_USFA ~ outcome_time,data = sum.data.ex.nons)
### Dunn test
dunnTest(DAG_USFA ~ outcome_time,
         data=sum.data.ex.nons,
         method="bh")    # Can adjust p-values;
#######statistical analysis For treatment-effect 

####kruskal.test
table(sum.data$CER_SFA)
sum.data.ex.nons<-sum.data[which(sum.data$outcome_time %in% c("Resolving_0","Non-resolving_0","Early-Nonsurvivors_0","Baseline_Baseline")),] %>% droplevels()

kruskal.test(LCER_USFA ~ TREATMENT,data = sum.data.ex.nons)
### Dunn test
dunnTest(LCER_USFA ~ TREATMENT,
         data=sum.data.ex.nons,
         method="bh")    # Can adjust p-values;
###
install.packages(pkgs="pastecs")
library(pastecs)

Idents(data_class)<-data_class$outcome
HC.class.data<-FetchData(subset(data_class,idents = c("Baseline")),vars = rownames(data_class) )
HC.class.data<-HC.class.data %>% mutate(CER_Hex= CER+HCER) %>% mutate(Total=CER+HCER+LCER+DCER+SM+PI+PE+PC+LPE+LPC+CE+MAG+TAG+DAG)
HC.class.data.stat<-stat.desc(HC.class.data)
HC.class.data$CER
data_class
