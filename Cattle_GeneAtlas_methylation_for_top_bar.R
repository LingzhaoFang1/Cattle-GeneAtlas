multimerge=function(mypath){
filenames=list.files(path=mypath,full.names=TRUE)
datalist=lapply(filenames,function(x)
{read.table(file=x,header=T)})
Reduce(function(x,y){merge (x,y,by="gene")},datalist)
}
library(tidyr)
name<-commandArgs(T)[1]
mydata=multimerge("//top10methylation//methylation")
mydata<-mydata[order(-mydata$value),]
rownames(mydata)<-mydata$gene
mydata<-mydata[,-1]
mydata$value<-NULL
a<-colMeans(mydata)
b<-sapply(mydata,function(x)sd(x)/sqrt(length(x)))
a<-as.data.frame(a)
b<-as.data.frame(b)
names<-as.data.frame(rownames(a))
colnames(names)<-"names"
rownames(names)<-names$names
tissue <-separate(names,names,into=c("tissue"),extra='drop')
bar<-as.data.frame(cbind(tissue,names,a,b))
colnames(bar)<-c("tissue","id","mean","se")

bar$color[bar$tissue=="Uterus"]<-"blue"
bar$color[bar$tissue=="Spleen"]<-"green"
bar$color[bar$tissue=="Rumen"]<-"orange"
bar$color[bar$tissue=="Ovary"]<-"yellow"
bar$color[bar$tissue=="Mammary"]<-"red"
bar$color[bar$tissue=="Lung"]<-"purple"
bar$color[bar$tissue=="Liver"]<-"pink"
bar$color[bar$tissue=="Muscle"]<-"violet"
bar$color[bar$tissue=="Kidney"]<-"black"
bar$color[bar$tissue=="Ileum"]<-"brown"
bar$color[bar$tissue=="Heart"]<-"gray"
bar$color[bar$tissue=="Frontal"]<-"#426F42"
bar$color[bar$tissue=="Lymphocyte"]<-"gold"
bar$color[bar$tissue=="Adipose"]<-"#FF1493"

bar<-bar[order(bar$mean,decreasing=TRUE),]
bar$id<-factor(bar$id,levels=bar$id)
library(ggplot2)
pdf(paste(name,'.top500.filtered.barplot.promoter-1500.pdf'),height=6,width=3.5)
ggplot(bar)+
geom_bar(aes(x=id,y=mean),fill=bar$color, color=bar$color,stat="identity",alpha=0.5)+xlab("")+ylab("Ajusted methylation level")+
geom_errorbar( aes(x=id,ymin=mean-se,ymax=mean+se),width=0.4,colour="black",alpha=0.9,size=0.4)+
coord_flip()+
 theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.off()
