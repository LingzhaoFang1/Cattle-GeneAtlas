#################################################################################################
#################################################################################################
###########################cattle-GeneAtlas project##############################################
##clean up R environment
rm(list = ls())

##loading required libraries
if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}

if (!require("reshape2")) {
  install.packages("reshape2", dependencies = TRUE)
  library(reshape2)
}

if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}

if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

if (!require("factoextra")) {
  install.packages("factoextra", dependencies = TRUE)##PCA plot
  library(factoextra)
}

if (!require("Rtsne")) {
  install.packages("Rtsne", dependencies = TRUE)
  library(Rtsne)
}

##############################################################################
##read into gene expression (FPKM) among 91 tissues and cell types in cattle, 
##which was computed uniformly using our RNA-seq pepiline. 
##To get access to expression data before paper publish: lingzhao.fang@igmm.ed.ac.uk

##set your current working directory, and put expression data in the dir
setwd("C:/Users/Lingzhao Fang/Desktop/US_project/Multiple_tissues/Cattle_expression/Final_RNA-seq/")


Gene_Expression <- get(load("./Gene_Expression_Altlas_723.Rdata"))
head(Gene_Expression[,1:5])

Smaple_Summary <- get(load("./Smaples_Summary.Rdata"))
head(Smaple_Summary)


####summary of tissue and cell types
Tissues_s <- Smaple_Summary$Tissue_s
table(Tissues_s)##smaple size per tissue/cell type
length(unique(Tissues_s))##number of unique tissue/cell types

###################################################################################################
##number of expressed genes across tissues and cell types##
##expressed genes with averaged FPKM>0 among samples in the tissue/cell type## 
Gene_Expression_T_F <- Gene_Expression>0 #convert to a T/F matrix
head(Gene_Expression_T_F[,1:10])


#####read into the 24616 ensemble genes based on Ensembl release 94 (UMD3.1.1)
UMD3_1 <- as.data.frame(fread("C:/Users/Lingzhao Fang/Desktop/US_project/Gestation Length GWAS/Genome_Annotation_5_18_2018/Cattle_ensemble_gene_5_18_2018.txt",
                              header = T,sep = "\t",stringsAsFactors = F))
head(UMD3_1)
##group genes based on gene type
table(UMD3_1$Gene.type)
##sort genes to the same order as in the expression data
UMD3_1_order <- UMD3_1[order(UMD3_1$Gene.stable.ID,decreasing = F),]
head(UMD3_1_order)
##check the order with the above expression atlas
identical(UMD3_1_order$Gene.stable.ID,colnames(Gene_Expression))



##############################################
##count expressed genes for each gene types

Expressed_genes <- array(NA,dim = c(dim(Gene_Expression_T_F)[1],10))
row.names(Expressed_genes) <- row.names(Gene_Expression_T_F)
colnames(Expressed_genes) <- names(table(UMD3_1_order$Gene.type))
head(Expressed_genes)


for(i in 1:dim(Gene_Expression_T_F)[1]){
  x <- table(UMD3_1_order$Gene.type[Gene_Expression_T_F[i,]])
  y <- names(x)
  z <- as.vector(table(UMD3_1_order$Gene.type[Gene_Expression_T_F[i,]]))
  Expressed_genes[i,y] <- z
}

head(Expressed_genes);class(Expressed_genes)
Expressed_genes_df <- as.data.frame(Expressed_genes)


Expressed_genes_df$Tissue <- Tissues_s
Expressed_genes_df$Tissue_Class <- Smaple_Summary$Tissue_visa

head(Expressed_genes_df)
Expressed_genes_df[is.na(Expressed_genes_df)] <- 0 ##missing values for gene counts as 0
Expressed_genes_df$Samples <- row.names(Expressed_genes_df)


###class samples by tissues; 
###averaged number of expressed genes among all samples per tissue 
Expressed_genes_class <- array(NA,dim = c(91,10))
row.names(Expressed_genes_class) <- unique(Tissues_s)
colnames(Expressed_genes_class) <- names(table(UMD3_1_order$Gene.type))
head(Expressed_genes_class)

for(i in 1:91){
  x <- Expressed_genes_df[Expressed_genes_df$Tissue==unique(Tissues_s)[i],1:10]
  y <- apply(x, 2, mean)
  Expressed_genes_class[unique(Tissues_s)[i],] <- y
}

head(Expressed_genes_class)

Expressed_genes_class_df <- as.data.frame(Expressed_genes_class)
Expressed_genes_class_df$Tissue <- row.names(Expressed_genes_class_df)
str(Expressed_genes_class_df)


####order tissues based on tissue classes

new_lables <- Smaple_Summary$Tissue_visa
#new_lables[Smaple_Summary$Tissue_s=="Milk_cells"] <- "Milk_cells"
head(new_lables)
Tissue_calss <- unstack(Tissues_s,Tissues_s~as.factor(new_lables))
str(Tissue_calss)

colors = rainbow(length(unique(new_lables)))
colors[4] <- "black"
#colors <- c('red','green','blue','orange','cyan','darkturquoise','goldenrod3','darkorchid3','mediumorchid1','lightseagreen','black','black','black','black','black')

names(colors) = unique(new_lables)
Tissue_calss <- Tissue_calss[names(colors)]

##unique tissues in each tissue class
Tissue_calss <- lapply(Tissue_calss, unique)

Tissue_calss_length <- lapply(Tissue_calss,length)
str(Tissue_calss_length)

###reorder the y axis names, that is the tissue names
Expressed_genes_class_df$Tissue <- factor(Expressed_genes_class_df$Tissue,levels = rev(unlist(Tissue_calss)))

####
Expressed_genes_class_df_melt <- melt(data = Expressed_genes_class_df,id.vars = "Tissue")
head(Expressed_genes_class_df_melt)

str(Expressed_genes_class_df_melt)

##plot the number of expressed genes among tissues and cell types
##change the path to your output dir
tiff(file = "C:/Users/Lingzhao Fang/Desktop/US_project/Multiple_tissues/manuscript/Figures/Run2/expressed_genes.tiff",
     res = 300, width = 2000, height = 2400,compression = "lzw")
p <- ggplot(data=Expressed_genes_class_df_melt, aes(y=value, x=Tissue, fill=factor(variable))) + 
  geom_bar(stat="identity", colour="black") +
  theme_bw() +ylab("Expressed Gene Counts")+theme(axis.text.y = element_text(hjust = 1, colour = rev(rep(colors,unlist(Tissue_calss_length)))))+
  coord_flip(ylim=c(1,20000))+
  labs(fill="Gene Type")
p       
dev.off()



#############################################################################################################
##############################cluster samples based on the gene expression using PCA and t-SNE
###################################################################################################################
############################################################################################################

#################### RPFK was log2 transformed################################

Gene_Expression_log<-log2(Gene_Expression+0.25)

head(Gene_Expression_log[,1:10])

###second methods for PCA
library(factoextra)


Gene_Expression_variance <- Gene_Expression_log[ , apply(Gene_Expression_log, 2, var) != 0]
head(Gene_Expression_variance[,1:10])

res.pca <- prcomp(Gene_Expression_variance, scale = TRUE,center =T) 

###group tissues
groups <- as.factor(Smaple_Summary$Tissue_visa)
unique(groups)

##col by tissue classes
Labels_factor <- as.factor(Smaple_Summary$Tissue_visa)

tiff(file = "Yourpath/PCA_Samples_log2.tiff",
     res = 300, width = 2800, height = 2000,compression = "lzw")
fviz_pca_ind(res.pca,
             geom = c("point"),
             col.ind = Labels_factor, # color by groups
             palette = colors,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)
dev.off()

######################################################
##########using t-SNE (non-linear) instead PCA(linear)

Train <- Gene_Expression_log
Texts <- Smaple_Summary$Tissue_s

new_lables <- Smaple_Summary$Tissue_visa
colors = rainbow(length(unique(new_lables)))
colors[4] <- "black"
names(colors) = unique(new_lables)

###clustering
tsne <- Rtsne(Train, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)

###ploting
##choose your own dir to store the fig
tiff(file = "Yourpath/tSNE.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
plot(tsne$Y, t='n',ylab="t-SNE-2",xlab = "t-SNE-1")
text(tsne$Y, labels=Texts, col=colors[new_lables],cex=0.5)

legend("topright", inset=c(0,0),horiz = F,  title="Tissue category",  # location of the legend on the heatmap plot
       legend = unique(new_lables), # category labels
       col = unique(colors[new_lables]),  # color key
       lty= 1,             # line style
       lwd = 5, cex = 0.6            # line width
)
dev.off()
