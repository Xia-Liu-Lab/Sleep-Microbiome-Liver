#20261111
#Microbial cinnamoylglycine mediates the effects of sleep quality on hepatic lipid metabolism
###------------------------------------------------------------------------------------------

microbiome<- read.csv()
combat_CMG<-read.csv()

median_value <- median(microbiome$SMHI, na.rm = TRUE)
microbiome$MRS_group <- ifelse(
  microbiome$SMHI <median_value, 0, 1
)

combat_CMG<-merge(microbiome,combat_CMG,by='row.names')

combat_CMG$log_CMG<-log10(combat_CMG$combat_CMG)

#sleep group
summary(glm(log_CMG~sleepgroup+age+gender+smoking+drinking+DD+MET+medication+hypertension+diabetes+obesity+EPI_eGFR,data=combat_CMG))
#MRS group
summary(glm(log_CMG~MRS_group+age+gender+smoking+drinking+DD+MET+medication+hypertension+diabetes+obesity+EPI_eGFR,data=combat_CMG))
#hepatic steatosis
summary(glm(log_CMG~FLI30+age+gender+smoking+drinking+DD+MET+medication+hypertension+diabetes+obesity+EPI_eGFR,data=combat_CMG))
#---Fibroscan
combat_fibroscan<-read.csv()
a<-glm(CAP~log_CMG+age+gender+smoking+drinking+DD+MET+medication+hypertension+diabetes+obesity+EPI_eGFR,data=combat_fibroscan)
summary(a)
a<-glm(FAST~log_CMG+age+gender+smoking+drinking+DD+MET+medication+hypertension+diabetes+obesity+EPI_eGFR,data=combat_fibroscan)
summary(a)

#-------------------------------RCS curve and threshold value------------------------------
library(rms)
library(segmented)
combat_CMG$FLI30<-as.factor(combat_CMG$FLI30)
dd<- datadist(combat_CMG)
options(datadist = 'dd')
#
fit<-ols(FLI~rcs(combat_CMG,3),data=combat_CMG)
fl<-Predict(fit,combat_CMG=seq(0,70,by=0.001))
anova(fit)
#
set.seed()
segmented_model<-segmented(lm(yhat~combat_CMG,data=fl))
plot.segmented(segmented_model)
segmented_model$psi


#------------------------------validation cohort from northern China--------------------
library(survival)

data<-read.csv()
a<-coxph(Surv(follow_up,incident_FLD)~factor(cut_group)+age+sex+smoke+PA+DM+HBP+obesity+EPI_eGFR,data=data)
summary(a)
#KM 
fit<-survfit(Surv(follow_up,incident_FLD)~cut_group,data=data)
fit
km_sum <- summary(fit, censored = TRUE)
cumhaz_source <- data.frame(
  group = if (is.null(km_sum$strata)) {
    "Overall"
  } else {
    sub("^group=", "", as.character(km_sum$strata))
  },
  time = km_sum$time,
  n_risk = km_sum$n.risk,
  n_event = km_sum$n.event,
  n_censor = km_sum$n.censor,
  survival = km_sum$surv,
  cumulative_hazard = -log(km_sum$surv),
  cumulative_hazard_lower_95CI = -log(km_sum$upper),
  cumulative_hazard_upper_95CI = -log(km_sum$lower)
)

#---------------------------validation from eastern China------------------------
library(ppcor)
liver_data <- read.csv()
liver_data$steatosis_group <- factor(liver_data$steatosis_group)

summary(glm(CMG ~ steatosis_group + Age + Gender + Smoking + Drinking + BMI,
             data = liver_data, family = gaussian()))

covariates <- liver_data[c("Age", "Gender", "Smoking", "Drinking", "BMI")]

ppcor::pcor.test(liver_data$CMG, liver_data$SAF_score,
                 covariates, method = "spearman")

ppcor::pcor.test(liver_data$CMG, liver_data$Inflammation_grade,
                 covariates, method = "spearman")

