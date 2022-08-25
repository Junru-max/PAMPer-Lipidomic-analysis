library(Seurat)
library(dplyr)
library(tidyr)
library(foreign)
library(survival)
library(survminer)
library(table1)

dim(pdata.trauma)
##Table 1
###For healthy Subjects
PAMPer.imp$AGE<- PAMPer.imp$AGE %>% as.numeric.factor()
median(PAMPer.imp$AGE[which(PAMPer.imp$outcome=="Baseline")])
IQR(PAMPer.imp$AGE[which(PAMPer.imp$outcome=="Baseline")])
table(PAMPer.imp$GENDER[which(PAMPer.imp$outcome=="Baseline")])
length(which(PAMPer.imp$GENDER[which(PAMPer.imp$outcome=="Baseline")]=="M"))/length(PAMPer.imp$GENDER[which(PAMPer.imp$outcome=="Baseline")])

table(PAMPer.imp$GENDER[which(PAMPer.imp$TIME.POINT %in% c("Baseline","0"))],PAMPer.imp$outcome[which(PAMPer.imp$TIME.POINT %in% c("Baseline","0"))] %>% droplevels) %>% chisq.test()
kruskal.test(AGE ~ outcome,data = PAMPer.imp@meta.data,subset=which(PAMPer.imp$TIME.POINT %in% c("Baseline","0"))
             )
?kruskal.test
###Table1 for patients in lipidomics dataset
pdata.trauma<-Pamper.pdata %>% filter(PAMPID %in% PAMPer.imp$`PAMPER ID NUMBER`)
pdata.trauma$TREATMENT<-PAMPer.imp$TREATMENT[match(pdata.trauma$PAMPID,PAMPer.imp$`PAMPER ID NUMBER`)] %>% droplevels()
pdata.trauma$outcome<-as.character(pdata.trauma$outcome)
pdata.trauma$traumatic_brain_injury<-pdata.trauma$traumatic_brain_injury %>% droplevels()
pdata.trauma$race<-as.character(pdata.trauma$race)
pdata.trauma$race[is.na(pdata.trauma$race)]<-"Non-white"
pdata.trauma$ed_coagulopathy[which(pdata.trauma$ed_coagulopathy==1)]<-"Yes"


pdata.trauma$outcome<-factor(pdata.trauma$outcome, 
                             levels=c("Resolving","Non-resolving","Early-Nonsurvivors","p-value"))


rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- pdata.trauma[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- kruskal.test(y ~ pdata.trauma$outcome)$p.value
    } else {
      p <- chisq.test(table(y, droplevels(pdata.trauma$outcome)))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "MEDIAN (IQR)"=sprintf("%s (&plusmn; %s)", MEDIAN, IQR)))
}

my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.0f %%)", FREQ, PCT))))
}



table1(~ age+gender+race+iss+ais_head+traumatic_brain_injury+initial_GCS+PH_sbp_70+vitals_hr+injury_type+TREATMENT+PH_time
       +PH_CPR+PH_intubation+PH_blood+PH_crystalloid+PH_prbc
       +transfusion_24h+prbc_24h+plasma_24h+platelets_24h+crystalloid_24h+vaso_24h+INR
       +ed_coagulopathy+ALI+NI+MOF+mech_vent_days+icu_los+hospital_los| outcome,
       data=pdata.trauma, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F,
       render.continuous=my.render.cont, render.categorical="FREQ (PCTnoNA%)",topclass="Rtable1-zebra",render.missing=NULL)
###

####Comorbility and Mediation
PAMPer_Com<-read.csv(file = "rawdata/PAMPer_com.csv") %>% filter(ï..stnum %in% pdata.trauma$stnum)
stnum_Com<-PAMPer_Com$ï..stnum
PAMPer_Med<-read.csv(file = "rawdata/PAMPer_mediation.csv") %>% filter(ï..stnum %in% pdata.trauma$stnum)
stnum_Med<-PAMPer_Med$ï..stnum
PAMPer_Com <- data.frame(lapply(PAMPer_Com, function(x){gsub(x,pattern="^3",replacement=NA)}),stringsAsFactors = T) %>% mutate(stnum=stnum_Com)
PAMPer_Med <- data.frame(lapply(PAMPer_Med, function(x){gsub(x,pattern="^3",replacement=NA)}),stringsAsFactors = T) %>% mutate(stnum=stnum_Med)

pdata.trauma<-pdata.trauma %>% left_join(PAMPer_Com,by="stnum") %>% left_join(PAMPer_Med,by="stnum")

table(pdata.trauma$stnum %in% PAMPer_Com$stnum)
table1(~ com_smoker+com_alcoholism+com_hypertension+com_diabetes+com_arrhythmia+com_history_mi+com_congestive_hf +com_copd +com_liver_disease+com_renal_dysfunction+com_coagulopathy+
         meds_antihypertensive+meds_statins+meds_betablockers+meds_corticosteriods+meds_anticoagulants+meds_nsaid+meds_asa+meds_other_antiplatelet+meds_oral_contraceptives| outcome,
       data=pdata.trauma, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F,
       render.continuous=my.render.cont, render.categorical="FREQ (PCTnoNA%)",topclass="Rtable1-zebra")
###




###survival analysis for treatment effect
####K-P CURVE
surv_object <- Surv(time = pdata.trauma$t_30d_mort_h,event = pdata.trauma$t_30d_censor)
fit1 <- survfit(surv_object ~ TREATMENT, data = pdata.trauma)
ggsurvplot(fit1, data = pdata.trauma,  risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),pval = TRUE,
           xscale=24,xlim=c(0,720),break.time.by=24,censor=T)

###COX REGRESSION
fit.coxph <- coxph(surv_object ~age+gender+iss+TREATMENT, 
                   data = pdata.trauma)
ggforest(fit.coxph,data = pdata.trauma)


####K-P CURVE
pdata.trauma.latephase<-pdata.trauma[which(pdata.trauma$outcome != "Early-Nonsurvivors"),] %>% drop.levels()
surv_object <- Surv(time = pdata.trauma.latephase$trans_icu_los,event = pdata.trauma.latephase$trans_icu_los_censor)
fit1 <- survfit(surv_object ~ TREATMENT, data = pdata.trauma.latephase)
dim(pdata.trauma.latephase)
length(pdata.trauma.latephase$TREATMENT)

ggsurvplot(fit1, data = pdata.trauma.latephase,  risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),pval = TRUE,
           xscale=1,xlim=c(0,30),break.time.by=4,censor=T)


###COX REGRESSION
fit.coxph <- coxph(surv_object ~age+gender+iss+TREATMENT, 
                   data = pdata.trauma.latephase)
ggforest(fit.coxph,data = pdata.trauma)
