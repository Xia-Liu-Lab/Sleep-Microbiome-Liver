#20261111
#Microbial cinnamoylglycine mediates the effects of sleep quality on hepatic lipid metabolism
###------------------------------------------------------------------------------------------
library(DESeq2)

#----------------------------Primary hepatocytes-------------------------------
#read
count_data <- read.csv()
sample_info <- read.csv()

# Match samples
count_data <- count_data[
  , rownames(sample_info),
  drop = FALSE
]

# Set reference group
sample_info$Group <- factor(
  sample_info$Group,
  levels = c("OA", "CMG")
)

# DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_info,
  design = ~ Group
)
#
dds <- DESeq(dds)
#
res <- results(dds,contrast = c("Group", "CMG", "OA"),alpha=0.05)
DEG_all <- as.data.frame(res)
#write.csv(DEG_all,file='')

#---------------------------Validation using human liver RNA-seq data (GSE135251)-------------------------
#read
count_data <- read.csv() 
sample_info <- read.csv() 

# Match samples
count_data <- count_data[
  , rownames(sample_info),
  drop = FALSE
]

# Set reference group
sample_info$Group <- factor(
  sample_info$Group,
  levels = c( "control","disease")
)

# DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_info,
  design = ~ Group
)
#
dds <- DESeq(dds)
#
res <- results(dds,contrast = c("Group","disease","control"),alpha=0.05)
DEG_all <- as.data.frame(res)

#write.csv(DEG_all,file='')


