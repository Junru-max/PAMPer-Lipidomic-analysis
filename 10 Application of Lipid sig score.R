library(Seurat)
library(pROC)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(caret)
library(glmnet)
library(dplyr)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(pROC)
library(scales)
library(survival)
library(survminer)
library(coxme)
library(rstatix)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

set.seed(11)
##Fig6 BCDEFG sFig6
##lineplot of Lip.sig.score for PAMPer

ggline(PAMPer.imp@meta.data, x = "TIME.POINT", y = "Lipid.sig.score4", 
       add = c("median_mad"),size = 1,point.size = 3,
       color = "outcome")+theme(aspect.ratio=1,legend.position = "bottom")
####K-P CURVE for PAMPer
pdata.trauma.latephase$Lipid.sig.score<-PAMPer.imp.72h$Lipid.sig.score[match(pdata.trauma.latephase$PAMPID,PAMPer.imp.72h$`PAMPER ID NUMBER`)]
pdata.trauma.latephase$Lipid.sig.score2<-PAMPer.imp.72h$Lipid.sig.score2[match(pdata.trauma.latephase$PAMPID,PAMPer.imp.72h$`PAMPER ID NUMBER`)]
pdata.trauma.latephase$Lipid.sig.score3<-PAMPer.imp.72h$Lipid.sig.score3[match(pdata.trauma.latephase$PAMPID,PAMPer.imp.72h$`PAMPER ID NUMBER`)]
pdata.trauma.latephase$Lipid.sig.score4<-PAMPer.imp.72h$Lipid.sig.score4[match(pdata.trauma.latephase$PAMPID,PAMPer.imp.72h$`PAMPER ID NUMBER`)]

pdata.trauma.latephase$clusters<-PAMPer.imp.72h$clusters[match(pdata.trauma.latephase$PAMPID,PAMPer.imp.72h$`PAMPER ID NUMBER`)]
pdata.trauma.latephase$clusters2<-PAMPer.imp.72h$clusters2[match(pdata.trauma.latephase$PAMPID,PAMPer.imp.72h$`PAMPER ID NUMBER`)]

surv_object <- Surv(time = pdata.trauma.latephase$trans_icu_los,event = pdata.trauma.latephase$trans_icu_los_censor)
fit1 <- survfit(surv_object ~ clusters2, data = pdata.trauma.latephase)
dim(pdata.trauma.latephase)
length(pdata.trauma.latephase$TREATMENT)

ggsurvplot(fit1, data = pdata.trauma.latephase,  risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),pval = TRUE,
           xscale=1,xlim=c(0,30),break.time.by=3,censor=T)

###COX REGRESSION for PAMPer
fit.coxph <- coxph(surv_object ~age+gender+iss+severe_head+TREATMENT+Lipid.sig.score2, 
                   data = pdata.trauma.latephase)
summary(fit.coxph)

ggforest(fit.coxph,data = pdata.trauma.latephase)
###Mixed Effects Cox
fit.coxme <- coxme(surv_object ~age+gender+iss+severe_head+TREATMENT+Lipid.sig.score2+ (1|SiteID), 
                   data = pdata.trauma.latephase)
summary(fit.coxme)
confint(fit.coxme)
ggforest(fit.coxme,data = pdata.trauma.latephase)
###Frailty model
fit.coxfra <- coxph(surv_object ~age+gender+iss+severe_head+TREATMENT+Lipid.sig.score+frailty(SiteID, distribution="gamma"), 
                   data = pdata.trauma.latephase)
summary(fit.coxfra)
class(fit.coxfra)<-"coxph"
fit.coxfra$coefficients
fit.coxph2<-fit.coxph
fit.coxph2$coefficients<-fit.coxfra$coefficients
fit.coxph2$coefficients[1]<- -0.0179

coef <- as.data.frame(tidy(fit.coxph2, conf.int = TRUE))

ggforest(fit.coxph2,data = pdata.trauma.latephase)
summary(fit.coxph2)
summary(fit.coxph)
summary(fit.coxfra)
exp(coef(fit.coxfra))
exp(coef(fit.coxph))
exp(confint(fit.coxfra))
exp(confint(fit.coxph))


###
##lineplot of Lip.sig.score for covid1
ggline(covid.imp@meta.data, x = "TIME", y = "Lipid.sig.score4", 
       add = c("median_mad"),size = 1,point.size = 3,
       color = "outcome")+theme(aspect.ratio=1,legend.position = "bottom")


###LR for covid-1
Idents(covid.imp)<-covid.imp$TIME
colnames(covid.imp@meta.data)
covid.reg<-FetchData(subset(covid.imp,ident=c("Before","<48h","D2-D5","D6_D14")),vars = colnames(covid.imp@meta.data)[c(20,19,24,28:35,41:43,37,45)])
covid.reg$outcome<-factor(covid.reg$outcome) %>% droplevels()
covid.reg$CRP.i..mg.L<-as.numeric.factor(covid.reg$CRP.i..mg.L)
model <- glm(outcome~., data = covid.reg,family = "binomial")
summary(model)
plot_model(model, sort.est = TRUE)
plot_model(model, transform = NULL,value.offset = .3,show.values = TRUE)+theme(aspect.ratio = 1)

###ROC for covid-1
par(pty="s")
hue_pal()(4)

pre.roc<-roc(covid.reg$outcome,covid.reg$Lipid.sig.score4,plot = T,legacy.axes=T,col="#F8766D",lwd=4,print.auc=T)
h0.roc<-roc(covid.reg$outcome,covid.reg$Lymphocyte.count..Ã.109.L,plot = T,legacy.axes=T,col="#7CAE00",lwd=4,print.auc=T,add=T)
h24.roc<-roc(covid.reg$outcome,covid.reg$CRP.i..mg.L,plot = T,legacy.axes=T,col="#00BFC4",lwd=4,print.auc=T,add=T)
legend("bottomright",legend=c("Lipid.sig.score","Lymphocyte","CRP"),col =hue_pal()(3) ,lwd=4)

roc(covid.reg$outcome,covid.reg$Lipid.sig.score2,plot = T,legacy.axes=T,col="#F8766D",lwd=4,print.auc=T)
roc(covid.reg$outcome,covid.reg$Lipid.sig.score3,plot = T,legacy.axes=T,col="#F8766D",lwd=4,print.auc=T)

##lineplot of Lip.sig.score for MM
MM.imp$outcome
ggline(MM.imp@meta.data, x = "TIME_POINT", y = "Lipid.sig.score2", 
       add = c("median_mad"),size = 1,point.size = 3,
       color = "outcome")+theme(aspect.ratio=1,legend.position = "bottom")

###sUrvival analysis for MM
####K-P CURVE
Idents(MM.imp)<-MM.imp$TIME_POINT
MM.imp.72h<-subset(MM.imp,idents="D2-D5")
MM.imp.72h<-ScaleData(MM.imp.72h,features = rownames(MM.imp.72h))
MM.imp.72h$Lipid.sig.score<-apply(MM.imp.72h@assays$RNA@scale.data[Lipid.sig,],2,mean)
MM.imp.72h$Lipid.sig.score2<-apply(MM.imp.72h@assays$RNA@scale.data[Lipid.sig2,],2,mean)
MM.imp.72h$Lipid.sig.score3<-MM.imp.72h@assays$RNA@scale.data[Lipid.sig3,]
MM.imp.72h$Lipid.sig.score4<-apply(MM.imp.72h@assays$RNA@scale.data[Lipid.sig4,],2,mean)

pdata.MM<-read.csv(file = "Rawdata/pdata.MM.csv")

pdata.MM$Lipid.sig.score<-MM.imp.72h$Lipid.sig.score[match(pdata.MM$Patient,MM.imp.72h$SUBJECT_ID)]
pdata.MM$Lipid.sig.score2<-MM.imp.72h$Lipid.sig.score2[match(pdata.MM$Patient,MM.imp.72h$SUBJECT_ID)]
pdata.MM$Lipid.sig.score3<-MM.imp.72h$Lipid.sig.score3[match(pdata.MM$Patient,MM.imp.72h$SUBJECT_ID)]
pdata.MM$Lipid.sig.score4<-MM.imp.72h$Lipid.sig.score4[match(pdata.MM$Patient,MM.imp.72h$SUBJECT_ID)]



pdata.MM<-pdata.MM[order(pdata.MM$Lipid.sig.score),]
pdata.MM$clusters<-c(rep("Low",29),rep("Med",29),rep("High",28))
pdata.MM<-pdata.MM[order(pdata.MM$Lipid.sig.score2),]
pdata.MM$clusters2<-c(rep("Low",29),rep("Med",29),rep("High",28))
table(pdata.MM$clusters,pdata.MM$clusters2)
surv_object <- Surv(time = pdata.MM$ICU.LOS)
fit1 <- survfit(surv_object ~ clusters2, data = pdata.MM)
ggsurvplot(fit1, data = pdata.MM,  risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),pval = TRUE,
           xscale=1,xlim=c(0,27),break.time.by=3,censor=F)
###COX REGRESSION
pdata.MM$GENDER<-factor(pdata.MM$GENDER,levels = c("Male","Female"))
pdata.MM$Head.Injury<-factor(pdata.MM$Head.Injury,levels = c(0,1),labels = c("No","Yes"))
pdata.MM$clusters<-factor(pdata.MM$clusters,levels = c("DA.Low","DA.Med","DA.High"))

fit.coxph <- coxph(surv_object ~Age+GENDER+ISS+Head.Injury+Lipid.sig.score3, 
                   data = pdata.MM)
ggforest(fit.coxph, data = pdata.MM)

##lineplot of Lip.sig.score for covid2

ggline(covid2.imp@meta.data, x = "GROUP.TIME", y = "Lipid.sig.score2", 
       add = c("median_mad"),size = 1,point.size = 3,
       color = "GROUP.TIME")+theme(aspect.ratio=1,legend.position = "bottom")
table(biomarker.heatmap$clusters,biomarker.heatmap$outcome2)


########statistical analysis for PAMPer
sum.data.ex.nons<-PAMPer.imp@meta.data[which(PAMPer.imp$outcome %in% c("Resolving","Non-resolving")),] %>% droplevels()
sum.data.ex.nons %>%
  group_by(outcome,TIME.POINT) %>%
  shapiro_test(Lipid.sig.score)
ggqqplot(sum.data.ex.nons, "Lipid.sig.score2", ggtheme = theme_bw()) +
  facet_grid(TIME.POINT ~ outcome, labeller = "label_both")
sum.data.ex.nons<-tibble(sum.data.ex.nons)
res.aov <-sum.data.ex.nons %>% anova_test(
  Lipid.sig.score2 ~ outcome*TIME.POINT )
#####2-way anova
get_anova_table(res.aov)

model <- lm(Lipid.sig.score2 ~ TIME.POINT * outcome, data = sum.data.ex.nons)
sum.data.ex.nons %>%
  group_by(TIME.POINT) %>%
  anova_test(Lipid.sig.score2 ~ outcome, error = model)
##
# pairwise comparisons
pwc <- sum.data.ex.nons %>% 
  group_by(TIME.POINT) %>%
  emmeans_test(Lipid.sig.score2 ~ outcome, p.adjust.method = "BH") 
pwc

####kruskal.test
sum.data.ex.nons<-PAMPer.imp@meta.data[which(PAMPer.imp$outcome_time %in% c("Resolving_0","Non-resolving_0","Early-Nonsurvivors_0")),] %>% droplevels()

kruskal.test(Lipid.sig.score2 ~ outcome_time,data = sum.data.ex.nons)
### Dunn test
dunnTest(Lipid.sig.score2 ~ outcome_time,
         data=sum.data.ex.nons,
         method="bh")    # Can adjust p-values;
###Group based heatmap

########statistical analysis for Covid-1

####kruskal.test
table(covid.imp$GROUP.TIME)

sum.data.ex.nons<-covid.imp@meta.data[which(covid.imp$GROUP.TIME %in% c("Non-severe_<48h","Non-severe_D2-D5","Non-severe_D6-D14","Severe_before_pro",
                                                                        "Severe_after_<48h","Severe_after_D6-D14")),] %>% droplevels()

kruskal.test(Lipid.sig.score2 ~ GROUP.TIME,data = sum.data.ex.nons)
### Dunn test
dunnTest(Lipid.sig.score2 ~ GROUP.TIME,
         data=sum.data.ex.nons,
         method="bh")    # Can adjust p-values;
###Group based heatmap

########statistical analysis for MM
table(MM.imp$TIME_POINT)
MM.imp$outcome
sum.data.ex.nons<-MM.imp@meta.data[which(MM.imp$outcome %in% c("Resolving","Non-resolving")),] %>% droplevels()
sum.data.ex.nons %>%
  group_by(outcome,TIME_POINT) %>%
  shapiro_test(Lipid.sig.score2)
ggqqplot(sum.data.ex.nons, "Lipid.sig.score2", ggtheme = theme_bw()) +
  facet_grid(TIME_POINT ~ outcome, labeller = "label_both")
sum.data.ex.nons<-tibble(sum.data.ex.nons)
res.aov <-sum.data.ex.nons %>% anova_test(
  Lipid.sig.score2 ~ outcome*TIME_POINT )
#####2-way anova
get_anova_table(res.aov)

model <- lm(Lipid.sig.score2 ~ TIME_POINT * outcome, data = sum.data.ex.nons)
sum.data.ex.nons %>%
  group_by(TIME_POINT) %>%
  anova_test(Lipid.sig.score2 ~ outcome, error = model)
##
# pairwise comparisons
pwc <- sum.data.ex.nons %>% 
  group_by(TIME_POINT) %>%
  emmeans_test(Lipid.sig.score2 ~ outcome, p.adjust.method = "BH") 
pwc

########statistical analysis for Covid-2

####kruskal.test
table(covid2.imp$GROUP.TIME)

sum.data.ex.nons<-covid2.imp@meta.data

kruskal.test(Lipid.sig.score2 ~ GROUP.TIME,data = sum.data.ex.nons)
### Dunn test
dunnTest(Lipid.sig.score2 ~ GROUP.TIME,
         data=sum.data.ex.nons,
         method="bh")    # Can adjust p-values;
###Group based heatmap
