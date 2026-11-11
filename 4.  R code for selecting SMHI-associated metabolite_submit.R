#20261111
#Microbial cinnamoylglycine mediates the effects of sleep quality on hepatic lipid metabolism
###------------------------------------------------------------------------------------------

# Packages loading
library(randomForest)
library(rfPermute)
library(mediation)

group_list<- read.csv()
metabolite<-read.csv()
metabolite_list<-metabolite[row.names(metabolite)%in%row.names(group_list),]

cor_data <- merge(group_list,metabolite_list,by="row.names")
row.names(cor_data) <- cor_data[,1]
cor_data <- cor_data[,-1]

# Stage 1 glm 
#GLM
results_diff <- data.frame(ID=NA,B=NA,P=NA,R2=NA,adjR2=NA)
n <- c(12:ncol(cor_data))
for (i in n) {
  model <- lm(cor_data[,i]~SMHI,data=cor_data)
  b <- summary(model)
  coef<-b$coefficients[2,1]
  P<-b$coefficients[2,4]
  R2<-b$r.squared
  adj.R2<-b$adj.r.squared
  diff_temp <- data.frame(ID=colnames(cor_data)[i],B=coef,P=P,R2=R2,adjR2=adj.R2)
  results_diff <- rbind(results_diff,diff_temp)
}
results_diff <- results_diff[-1,]
#write.csv(results_diff,...)

#Permutation
rho_result_permutated <- as.data.frame(array(dim = c(ncol(metabolite_list),1001)))
colnames(rho_result_permutated)[1] <- "ID"
rho_result_permutated$ID <- colnames(metabolite_list)
n <- c(1:ncol(metabolite_list))
p <- c(1:1000)
set.seed()
for (i in n) {
  cor_data_temp <- cor_data[which(cor_data[,i+11]>-1000),]
  for (j in p) {
    cor_data_temp_permutate <- cor_data_temp
    cor_data_temp_permutate[,i+11] <- sample(cor_data_temp_permutate[,i+11])
    model <- lm(cor_data_temp_permutate[,i+11]~SMHI,data=cor_data_temp_permutate)
    b <- summary(model)
    coef<-b$coefficients[2,1]
    P<-b$coefficients[2,4]
    rho_result_permutated[i,j+1] <-P
  }
}
#write.csv(rho_result_permutated,...)

rho_result <- results_diff

rho_result$FDR <- NA

row.names(rho_result)<-rho_result[,1]
row.names(rho_result_permutated)<-rho_result_permutated[,1]
rho_result<-rho_result[,-1]
rho_result_permutated<-rho_result_permutated[,-1]

n <- c(1:nrow(rho_result))
for (i in n) {
  fdr <- length(which(abs(rho_result_permutated[i,])<abs(rho_result$P[i])))/ncol(rho_result_permutated)
  rho_result$FDR[i] <- fdr
}

list<-subset(rho_result,FDR<0.05)

metabolite_table<- metabolite_list[,colnames(metabolite_list)%in%row.names(list)]
metabolite_table<-merge(group_list,metabolite_table,by="row.names")
#write.csv(list,...)
#write.csv(metabolite_table, ...)

#Stage 2 Random forest
metabolite_name<-row.names(list)
rf_table<-metabolite_table[,c("FLI30",metabolite_name)]
rf_table$FLI30<-as.factor(rf_table$FLI30)
set.seed()
rf<-randomForest(FLI30~.,data=rf_table)
rfP <- rfPermute(FLI30~., data =rf_table)
VIP<-as.data.frame(importance(rfP))
write.csv(VIP,...)

# Stage 3 glm 
results_diff <- data.frame(ID=NA,P=NA,B=NA,low=NA,high=NA)
n <- match(rownames(list), colnames(metabolite_table))
for (i in n) {
  model <- glm(FLI~ metabolite_table[,i]+age+gender+smoking+drinking+DD+MET+diabetes+hypertension+obesity+medication,data=metabolite_table)
  b <- summary(model)
  p <- b$coefficients[2,4]
  B <-b$coefficients[2,1]
  c<-confint(model)
  low<-c[2,1]
  high<-c[2,2]
  diff_temp <- data.frame(ID=colnames(metabolite_table)[i],P=p,B=B,low=low,high=high)
  results_diff <- rbind(results_diff,diff_temp)
}
results_diff <- results_diff[-1,]
#write(results_diff,...)



#------------------------------Correlation of pathways and metabolites----------------------------
pathway_data<-read.csv()
Clinic<-read.csv()
metabolite_inclusion<-read.csv()
metabolite_data<-metabolite_list[,rownames(metabolite_inclusion)]

##
harmonized<-merge(metabolite_data,pathway_data,by = "row.names")
harmonized$Row.names <- NULL
####cor.test for pathway and metabolite data
n <- c(49:56)
p <- c(1:48)
pval <- as.data.frame(array(dim = c(8,48)))
colnames(pval) <- colnames(harmonized)[c(1:48)]
row.names(pval) <- colnames(harmonized)[c(49:56)]
rval <- pval
for (j in p) {
  for (i in n) {
    a <- cor.test(harmonized[,j],harmonized[,i],method = "spearman")
    pval[colnames(harmonized)[i],colnames(harmonized)[j]]<-a$p.value
    rval[colnames(harmonized)[i],colnames(harmonized)[j]] <- a$estimate
  }
}
#write.csv(rval,)

#Permutation
pathway<- harmonized[,49:56]
metabolite<- harmonized[,c(1:48)]
#
new_data <- cbind(pathway,metabolite)
#
permutated_p <- as.data.frame(array(dim = c(ncol(pathway),ncol(metabolite))))
#
row.names(permutated_p) <- colnames(pathway)
colnames(permutated_p) <- colnames(metabolite)
n <- c(1:ncol(pathway))
p <- c((ncol(pathway)+1):ncol(new_data))
#Permutation
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
#write.csv(permutated_p,....)

#------------------------------Correlation of 12 metabolites and clinical indicators--------------
metabolite_inclusion<-read.csv()
metabolite_data<-metabolite_list[,rownames(metabolite_inclusion)]

##
harmonized<-merge(metabolite_data,Clinic,by = "row.names")
harmonized<-harmonized[,-c(22)]
harmonized$Row.names <- NULL
harmonized<-harmonized[,c(1:18,22:25,19:21)]
#
n <- c(13:25)
p <- c(1:12)
pval <- as.data.frame(array(dim = c(13,12)))
colnames(pval) <- colnames(harmonized)[c(1:12)]
row.names(pval) <- colnames(harmonized)[c(13:25)]
rval <- pval
for (j in p) {
  for (i in n) {
    a <- cor.test(harmonized[,j],harmonized[,i],method = "spearman")
    pval[colnames(harmonized)[i],colnames(harmonized)[j]]<-a$p.value
    rval[colnames(harmonized)[i],colnames(harmonized)[j]] <- a$estimate
  }
}
#write.csv(rval,...)

#Permutation
#
metadata<- harmonized[,13:25]
metabolite<- harmonized[,c(1:12)]
#
new_data <- cbind(metadata,metabolite)
#
permutated_p <- as.data.frame(array(dim = c(ncol(metadata),ncol(metabolite))))
#
row.names(permutated_p) <- colnames(metadata)
colnames(permutated_p) <- colnames(metabolite)
n <- c(1:ncol(metadata))
p <- c((ncol(metadata)+1):ncol(new_data))
#
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

#--------------------------mediation of metabolites-----------------------------------------------
microbiome<- read.csv()
metabolite_inclusion<-read.csv()
metabolite_data<-metabolite_list[,rownames(metabolite_inclusion)]
med<-merge(microbiome,metabolite_data,by="row.names")
row.names(med)<-med$Row.names
med<-med[,-c(1)]
#mediation
#Putting each potential metabolite into the model as a mediator
#for example CMG (cinnamoylglycine)
#continuous FLI
med.xm <- glm(CMG ~ SMHI + age+gender+smoking+drinking+DD+MET+hypertension+diabetes+obesity+medication, data = med)
med.xy <- glm(FLI~ SMHI +CMG + age+gender+smoking+drinking+DD+MET+hypertension+diabetes+obesity+ medication, data = med)
set.seed()
mod.med <- mediate(med.xm, med.xy, treat = 'SMHI', mediator = 'CMG', sims = 1000, boot = T)
summary(mod.med)
plot(mod.med)
#prevalent FLI
med$FLI30<-as.factor(med$FLI30)
med.xm <- glm(CMG ~ SMHI + age+gender+smoking+drinking+DD+MET+hypertension+diabetes+obesity+medication, data = med)
med.xy <- glm(FLI30~ SMHI +CMG + age+gender+smoking+drinking+DD+MET+hypertension+diabetes+obesity+ medication, family=binomial,data = med)
mod.med <- mediate(med.xm, med.xy, treat = 'SMHI', mediator = 'CMG', sims = 1000, boot = T)
summary(mod.med)
plot(mod.med)
