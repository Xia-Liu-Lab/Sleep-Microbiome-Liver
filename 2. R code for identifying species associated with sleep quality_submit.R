#20261111
#Microbial cinnamoylglycine mediates the effects of sleep quality on hepatic lipid metabolism
###------------------------------------------------------------------------------------------

# Packages loading
options("scipen"=999, "digits"=4)
library(glmnet)
library(mediation)

#### read
group<-read.csv() #metadata
species_raw_table<-read.csv() #species file

#### Prevalence filtering
delete_list <- c()
n <- c(1:ncol(species_raw_table))
for (i in n) {
  if(length(which(species_raw_table[,i]>0))<(nrow(species_raw_table)*0.01)){
    delete_list <- append(delete_list,i)
  }  
}
species_prev_filter <- species_raw_table[,-delete_list]

# Stage 1 Spearman correlation
cor_data <- merge(group,species_prev_filter,by="row.names")
row.names(cor_data) <- cor_data[,1]
cor_data <- cor_data[,-1]
#
rho_result <- as.data.frame(array(dim = c(ncol(species_prev_filter),3)))
colnames(rho_result) <- c("ID","rho","p")
rho_result$ID <- colnames(species_prev_filter)
n <- c(1:ncol(species_prev_filter))
for (i in n) {
  cor_data_temp <- cor_data[which(cor_data[,i+17]>0),]
  a <- cor.test(cor_data_temp$sleepscore,cor_data_temp[,i+17],method = "spearman")
  rho <- a$estimate
  p <- a$p.value
  rho_result$rho[i] <- rho
  rho_result$p[i] <- p
}
#Permutation
rho_result_permutated <- as.data.frame(array(dim = c(ncol(species_prev_filter),1001)))
colnames(rho_result_permutated)[1] <- "ID"
rho_result_permutated$ID <- colnames(species_prev_filter)
n <- c(1:ncol(species_prev_filter))
p <- c(1:1000)
set.seed()
for (i in n) {
  cor_data_temp <- cor_data[which(cor_data[,i+17]>0),]
  for (j in p) {
    cor_data_temp_permutate <- cor_data_temp
    cor_data_temp_permutate[,i+17] <- sample(cor_data_temp_permutate[,i+17])
    a <- cor.test(cor_data_temp_permutate$sleepscore,cor_data_temp_permutate[,i+17],method = "spearman")
    rho <- a$p.value
    rho_result_permutated[i,j+1] <- rho
  }
}
#calculate FDR
rho_result$FDR <- NA
row.names(rho_result)<-rho_result[,1]
row.names(rho_result_permutated)<-rho_result_permutated[,1]
rho_result<-rho_result[,-1]
rho_result_permutated<-rho_result_permutated[,-1]

n <- c(1:nrow(rho_result))
for (i in n) {
  fdr <- length(which(abs(rho_result_permutated[i,])<abs(rho_result$p[i])))/ncol(rho_result_permutated)
  rho_result$FDR[i] <- fdr
}
#Selecting species with FDR < 0.1
rho_result <- rho_result[which(rho_result$FDR<0.1),]
specie_cor <- species_prev_filter[,colnames(species_prev_filter)%in%row.names(rho_result)]
specie_cor <- merge(group,specie_cor,by="row.names")
row.names(specie_cor) <- specie_cor$Row.names
specie_cor <- specie_cor[,-grep("Row",colnames(specie_cor))]

# Stage 2 LASSO
lasso_model <- cv.glmnet(x=as.matrix(specie_cor[18:ncol(specie_cor)]), y=specie_cor$sleepscore, alpha=1, family = 'gaussian', type.measure='deviance',nfolds = 5)
Coefficients_temp <- as.matrix(coef(lasso_model,s=lasso_model$lambda.min))
Coefficients_temp <- data.frame(Coefficients_temp)
colnames(Coefficients_temp)[1] <- paste("round",1,sep = "")
lasso_permutated <- Coefficients_temp
#
n <- c(2:1000)
set.seed()
for (i in n) {
  lasso_model <- cv.glmnet(x=as.matrix(specie_cor[18:ncol(specie_cor)]), y=specie_cor$sleepscore, alpha=1, family = 'gaussian', type.measure='deviance',nfolds = 5)
  Coefficients_temp <- as.matrix(coef(lasso_model,s=lasso_model$lambda.min))
  Coefficients_temp <- data.frame(Coefficients_temp)
  colnames(Coefficients_temp)[1] <- paste("round",i,sep = "")
  lasso_permutated <- cbind(lasso_permutated,Coefficients_temp)
}
lasso_permutated <- lasso_permutated[-1,] 
lasso_select <- as.data.frame(array(dim = c(nrow(lasso_permutated),2)))
colnames(lasso_select) <- c("ID","frequency")
lasso_select$ID <- row.names(lasso_permutated)
n <- c(1:nrow(lasso_permutated))
for (i in n) {
  freq <- length(which(lasso_permutated[i,]!=0))/ncol(lasso_permutated)  
  lasso_select$frequency[i] <- freq
}
#Selecting Speices with Frequency >50% 
lasso_list <- lasso_select$ID[which(lasso_select$frequency>0.5)]
lasso_table <- species_prev_filter[,colnames(species_prev_filter)%in%lasso_list]
lasso_table <- merge(group,lasso_table,by="row.names")
row.names(lasso_table) <- lasso_table$Row.names
lasso_table <- lasso_table[,-grep("Row",colnames(lasso_table))]

#-------Variance in sleep quality explained by each microbial species --------------
results_diff <- data.frame(ID=NA,P=NA,R2=NA,adj.R2=NA)
n <- c(18:31)
for (i in n) {
  model <- lm(sleepscore~lasso_table[,i]+age+gender+BMI,data=lasso_table)
  b <- summary(model)
  b
  R2<-b$r.squared
  adj.R2<-b$adj.r.squared
  fvalue <- b$fstatistic
  pvalue<- 1-pf(fvalue[1],fvalue[2],fvalue[3])
  diff_temp <- data.frame(ID=colnames(lasso_table)[i],P=pvalue,R2=R2,adj.R2=adj.R2)
  results_diff <- rbind(results_diff,diff_temp)
}
results_diff <- results_diff[-1,]


# ----------Construction of SMHI--------------------
#Step 1 calculate weights for the assocaition of sleep quality
species_coeff <- data.frame(ID=NA,B=NA,low=NA,high=NA)
n <- c(1:(ncol(lasso_table)-17))
for (i in n) {
  model <- glm(sleepscore~lasso_table[,i+17],data = lasso_table)
  b <- summary(model)
  B <-b$coefficients[2,1]
  c<-confint(model)
  low<-c[2,1]
  high<-c[2,2]
  diff_temp <- data.frame(ID=colnames(lasso_table)[i+17],B=B,low=low,high=high)
  species_coeff<- rbind(species_coeff,diff_temp)
}
species_coeff <- species_coeff[-1,]
#Step 2 construction of SMHI
lasso_table$SMHI <- NA
n <- c(1:nrow(lasso_table))
for (i in n) {
  lasso_table$SMHI[i] <- lasso_table[i,18]*species_coeff$B[1]+
    lasso_table[i,19]*species_coeff$B[2]+
    lasso_table[i,20]*species_coeff$B[3]+
    lasso_table[i,21]*species_coeff$B[4]+
    lasso_table[i,22]*species_coeff$B[5]+
    lasso_table[i,23]*species_coeff$B[6]+
    lasso_table[i,24]*species_coeff$B[7]+
    lasso_table[i,25]*species_coeff$B[8]+
    lasso_table[i,26]*species_coeff$B[9]+
    lasso_table[i,27]*species_coeff$B[10]+
    lasso_table[i,28]*species_coeff$B[11]+
    lasso_table[i,29]*species_coeff$B[12]+
    lasso_table[i,30]*species_coeff$B[13]+
    lasso_table[i,31]*species_coeff$B[14]
}

#-------- Association with hepatic steatosis----
correlation<-lasso_table
correlation$FLI30<-as.factor(correlation$FLI30)
correlation$SMHI_quartile <- cut(
  correlation$SMHI,
  breaks = quantile(
    correlation$SMHI,
    probs = seq(0, 1, 0.25),
    na.rm = TRUE
  ),
  include.lowest = TRUE,
  labels = c("Q1", "Q2", "Q3", "Q4")
)
table(correlation$SMHI_quartile)

summary(glm(FLI30~ factor(SMHI_quartile)+age+gender+smoking+drinking+DD+MET+diabetes+hypertension+obesity+medication,family=binomial(),data=correlation))
confint(glm(FLI30~ factor(SMHI_quartile)+age+gender+smoking+drinking+DD+MET+diabetes+hypertension+obesity+medication,family=binomial(),data=correlation))

summary(glm(FLI~ factor(SMHI_quartile)+age+gender+smoking+drinking+DD+MET+diabetes+hypertension+obesity+medication,data=correlation))
confint(glm(FLI~ factor(SMHI_quartile)+age+gender+smoking+drinking+DD+MET+diabetes+hypertension+obesity+medication,data=correlation))

#------------ mediation ----------
media <- lasso_table

#continuous outcome
med.xm <- glm(SMHI ~ sleepscore + age+gender+smoking+drinking+DD+MET+medication+hypertension+diabetes+obesity, data = media)
med.xy <- glm(FLI~ sleepscore + SMHI + age+gender+smoking+drinking+DD+MET+medication+hypertension+diabetes+obesity, data =media)
set.seed()
mod.med <- mediate(med.xm, med.xy, treat = 'sleepscore', mediator = 'SMHI', sims = 1000, boot = T)
summary(mod.med)

# hepatic steatosis
med.xm <- glm(SMHI ~ sleepscore + age+gender+smoking+drinking+DD+MET+medication+hypertension+diabetes+obesity, data = media)
med.xy <- glm(FLI30~ sleepscore + SMHI + age+gender+smoking+drinking+DD+MET+medication+hypertension+diabetes+obesity, family=binomial(),data = media)
set.seed()
mod.med <- mediate(med.xm, med.xy, treat = 'sleepscore', mediator = 'SMHI', sims = 1000, boot = T)
summary(mod.med)

#---------------Correlation of species with clinic indicators----------------
#read 
microbiome<- read.csv()
Clinic<-read.csv()
harmonized<-merge(microbiome,Clinic,by.x = "row.names",by.y = "Row.names")
harmonized$Row.names <- NULL
#cor.test for microbiome and clinic data
#
n <- c(1:15)
p <- c(16:21,26:29,22,23,25,31,30)
pval <- as.data.frame(array(dim = c(15,15)))
colnames(pval) <- colnames(harmonized)[c(16:21,26:29,22,23,25,31,30)]
row.names(pval) <- colnames(harmonized)[1:15]
rval <- pval
for (j in p) {
  for (i in n) {
    a <- cor.test(harmonized[,j],harmonized[,i],method = "spearman")
    pval[colnames(harmonized)[i],colnames(harmonized)[j]]=a$p.value
    rval[colnames(harmonized)[i],colnames(harmonized)[j]] <- a$estimate
  }
}
#write.csv(rval,...)
# permutation FDR
microbiome<- harmonized[,1:15]
metadata<- harmonized[,c(16:21,26:29,22,23,25,31,30)]
new_data <- cbind(microbiome,metadata)
#
permutated_p <- as.data.frame(array(dim = c(ncol(microbiome),ncol(metadata))))
#
row.names(permutated_p) <- colnames(microbiome)
colnames(permutated_p) <- colnames(metadata)
n <- c(1:ncol(microbiome))
p <- c((ncol(microbiome)+1):ncol(new_data))
#permutation
set.seed()
round <- c(1:1000)
for (i in n) {
  for (j in p) {
    raw_p <- as.numeric(cor.test(new_data[,i],new_data[,j],method = "spearman")$p.value)
    list <- c()
    for (k in round) {
      temp <- new_data
      temp[,i] <- sample(temp[,i])
      p_val <- as.numeric(cor.test(temp[,i],temp[,j],method = "spearman")$p.value)
      list <- append(list,p_val)
    }
    permutated_value <- length(which(list<raw_p))/length(round)
    permutated_p[i,(j-length(n))] <- permutated_value
  }
}
#write.csv(permutated_p,...)