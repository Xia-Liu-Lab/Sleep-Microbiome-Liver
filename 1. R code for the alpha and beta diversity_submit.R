#20261111
#Microbial cinnamoylglycine mediates the effects of sleep quality on hepatic lipid metabolism
###------------------------------------------------------------------------------------------

library(ggplot2)
library(vegan)
library(gghalves)
library(ape)
#---alpha diversity---
# file read & reshape
pri <- read.csv()  # Genus file
pri <- as.data.frame(t(pri))
# shannon calculate
shannon <- data.frame(diversity(pri,"shannon"))
names(shannon) <- "shannon"
shannon$MetaID <- row.names(shannon)
#add group
group <- read.csv() #metadata file
harmo <- merge(shannon, group, by.x = "MetaID",by.y="row.names")
harmo$sleepgroup<-as.factor(harmo$sleepgroup)
#plot
p<-ggplot(harmo,aes(x= sleepgroup, y= shannon,fill=sleepgroup,color=sleepgroup))+
  scale_fill_manual(values=c('#C2202F','#E3AB5A','#3172A7'))+
  scale_color_manual(values=c('#C2202F','#E3AB5A','#3172A7'))
p
#
p1<-p+geom_half_violin(position=position_nudge(x=0.1,y=0),
                       side='R',adjust=1.2,trim=F,color=NA,alpha=0.6)
p1
#
p2<-p1+geom_half_point(position=position_nudge(x=-0.35,y=0),size =0.5, 
                       shape =19,range_scale = 0.5,alpha=0.3)
p2
#
p3<-p2+geom_boxplot(outlier.shape = NA, #隐藏离群点；
                    width =0.1,
                    alpha=0.5)
p3
#
p4<-p3+coord_flip()
p4
#
p5 <- p4+theme_bw()+theme(panel.grid=element_blank(),legend.position = 'none',
                          text = element_text(size = 8,family='sans'),
                          axis.text = element_text(size = 8,family='sans',color='#000000'),
                          axis.title = element_text(size = 8,family='sans',color='#000000'))+
  labs(y='Shannon Index',x='')
p5

#---beta diversity---
#read
pri <- read.csv()  # Genus file
pri <- as.data.frame(t(pri))
group <- read.csv() #group file
pri<-pri[row.names(pri)%in%row.names(group),]
group$sleepgroup<-as.factor(group$sleepgroup)
group <- group$sleepgroup

#PCOA
dist_bray <- vegdist(pri, method = "bray")
pcoa_result <- pcoa(dist_bray)
scores <- data.frame(
  Sample = rownames(pri),
  PCo1 = pcoa_result$vectors[, 1],
  PCo2 = pcoa_result$vectors[, 2],
  Group = group
)

#plot
ggplot(scores, aes(x = PCo1, y = PCo2, color = factor(Group))) +
  geom_point(size = 1) +
  theme_bw()+
  theme(panel.grid=element_blank(),legend.position = 'none',
        text = element_text(size = 8,family='sans'),
        axis.text = element_text(size = 8,family='sans',color='#000000'),
        axis.title = element_text(size = 8,family='sans',color='#000000'))+
  scale_color_manual(values = c('#D96D7B','#E7BD7C','#5191B6'))+
  stat_ellipse(level = 0.95,linetype=2) +
  labs(x = paste("PCo1 (", round(pcoa_result$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste("PCo2 (", round(pcoa_result$values$Relative_eig[2] * 100, 2), "%)"))


