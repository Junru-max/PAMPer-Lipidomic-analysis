library(Seurat)
library(dplyr)
library(tidyr)
library(foreign)
library(table1)

##sTable1
######Table for LRS groups
pdata.trauma.latephase<-pdata.trauma %>% filter(PAMPID %in% PAMPer.imp.72h$`PAMPER ID NUMBER`)
pdata.trauma.latephase$clusters<-PAMPer.imp.72h$clusters[match(pdata.trauma.latephase$PAID,PAMPer.imp.72h$`PAMPER ID NUMBER`)]
pdata.trauma.latephase$clusters<-factor(pdata.trauma.latephase$clusters, 
                               levels=c("Low","Med","High","P-value"), 
                               labels=c("LRS-Low","LRS-Med","LRS-High","P-value"))
pdata.trauma.latephase$outcome<-pdata.trauma.latephase$outcome %>% droplevels()
pdata.trauma.latephase$traumatic_brain_injury<-pdata.trauma.latephase$traumatic_brain_injury %>% droplevels()
pdata.trauma.latephase$race<-as.character(pdata.trauma.latephase$race)
pdata.trauma.latephase$race[pdata.trauma.latephase$race!="1:White"]<-"Non-white"
pdata.trauma.latephase$ed_coagulopathy[which(pdata.trauma.latephase$ed_coagulopathy==1)]<-"Yes"


rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- pdata.trauma.latephase[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- kruskal.test(y ~ pdata.trauma.latephase$clusters)$p.value
    } else {
      p <- chisq.test(table(y, droplevels(pdata.trauma.latephase$clusters)))$p.value
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


table1(~  age+gender+race+iss+ais_head+traumatic_brain_injury+initial_GCS+PH_sbp_70+vitals_hr+injury_type+TREATMENT+PH_time
       +PH_CPR+PH_intubation+PH_blood+PH_crystalloid+PH_prbc
       +transfusion_24h+prbc_24h+plasma_24h+platelets_24h+crystalloid_24h+vaso_24h+INR
       +ed_coagulopathy+ALI+NI+MOF+mech_vent_days+icu_los+hospital_los+mortality_30d+outcome| clusters,
       data=pdata.trauma.latephase, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F,
       render.continuous=my.render.cont, render.categorical="FREQ (PCTnoNA%)",topclass="Rtable1-zebra",render.missing=NULL)


