#20261111
#Microbial cinnamoylglycine mediates the effects of sleep quality on hepatic lipid metabolism
###------------------------------------------------------------------------------------------

# Packages loading
library(randomForest)
library(rfPermute)
#read 
microbiome<- read.csv()
pathway_raw_table <- read.csv()

#
delete_list <- c()
n <- c(1:ncol(pathway_raw_table))
for (i in n) {
  if(length(which(pathway_raw_table[,i]>0))<(nrow(pathway_raw_table)*0.01)){
    delete_list <- append(delete_list,i)
  }  
}
pathway_prev_filter <- pathway_raw_table[,-delete_list]

# Stage 1 glm 
group_list<-read.csv()
cor_data <- merge(group_list,pathway_prev_filter,by="row.names")
row.names(cor_data) <- cor_data[,1]
cor_data <- cor_data[,-1]

results_diff <- data.frame(ID=NA,B=NA,P=NA,R2=NA,adjR2=NA)
n <- 12:ncol(cor_data)
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
#Permutation
rho_result_permutated <- as.data.frame(array(dim = c(ncol(pathway_prev_filter),1001)))
colnames(rho_result_permutated)[1] <- "ID"
rho_result_permutated$ID <- colnames(pathway_prev_filter)
n <- c(1:ncol(pathway_prev_filter))
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
rho_result_permutated

rho_result<-results_diff
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
pathway_list<-subset(rho_result,FDR<0.05)
#write.csv(pathway_list,...)

#Stage 2 Random forest 
rf_table <- read.csv()  #only FLI30 and potential pathway
rf_table$FLI30<-as.factor(rf_table$FLI30)
set.seed()
rf<-randomForest(FLI30~.,data=rf_table)
rfP <- rfPermute(FLI30~., data =rf_table)
VIP<-as.data.frame(importance(rfP))
#write.csv(VIP,....)

#Stage 3 glm 
pathway_table <- read.csv() #metadata and potential pathway
results_diff <- data.frame(ID=NA,P=NA,B=NA,low=NA,high=NA)
n <- c(48:ncol(pathway_table))
for (i in n) {
  model <- glm(FLI~ pathway_table[,i]+age+gender+smoking+drinking+DD+MET+diabetes+hypertension+obesity+medication,data=pathway_table)
  b <- summary(model)
  p <- b$coefficients[2,4]
  B <-b$coefficients[2,1]
  c<-confint(model)
  low<-c[2,1]
  high<-c[2,2]
  diff_temp <- data.frame(ID=colnames(pathway_table)[i],P=p,B=B,low=low,high=high)
  results_diff <- rbind(results_diff,diff_temp)
}
results_diff <- results_diff[-1,]
#write.csv(results_diff,...)

#------ Correlation of pathway and clinical indicators-------
pathway<- read.csv()
pathway<-pathway[,c(19:26)]
Clinic<-read.csv()
harmonized<-merge(Clinic,pathway,by="row.names")
harmonized$Row.names <- NULL

#cor.test for pathway and clinic data
n <- c(15:22)
p <- c(1:6,11:14,7:8,10)
pval <- as.data.frame(array(dim = c(8,13)))
colnames(pval) <- colnames(harmonized)[c(1:6,11:14,7:8,10)]
row.names(pval) <- colnames(harmonized)[15:22]
rval <- pval
for (j in p) {
  for (i in n) {
    a <- cor.test(harmonized[,j],harmonized[,i],method = "spearman")
    pval[colnames(harmonized)[i],colnames(harmonized)[j]]=a$p.value
    rval[colnames(harmonized)[i],colnames(harmonized)[j]] <- a$estimate
  }
}
#write.csv(rval,...)
#permutation
pathway<- harmonized[,15:22]
metadata<- harmonized[,c(1:6,11:14,7:8,10)]
#
new_data <- cbind(pathway,metadata)
#
permutated_p <- as.data.frame(array(dim = c(ncol(pathway),ncol(metadata))))
#
row.names(permutated_p) <- colnames(pathway)
colnames(permutated_p) <- colnames(metadata)
n <- c(1:ncol(pathway))
p <- c((ncol(pathway)+1):ncol(new_data))
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


