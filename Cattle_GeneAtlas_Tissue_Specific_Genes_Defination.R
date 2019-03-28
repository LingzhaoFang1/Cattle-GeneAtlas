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

if (!require("stats")) {
  install.packages("stats", dependencies = TRUE)
  library(stats)
}

##############################################################################
##read into gene expression (FPKM) among 91 tissues and cell types in cattle, 
##which was computed uniformly using our RNA-seq pepiline. 
##To get access to expression data before paper publish: lingzhao.fang@igmm.ed.ac.uk

##set your current working directory, and put expression data in the dir
setwd("C:/Users/Lingzhao Fang/Desktop/US_project/Multiple_tissues/Cattle_expression/Final_RNA-seq/")


Gene_Expression <- get(load("./Gene_Expression_Altlas_723.Rdata"))
head(Gene_Expression[,1:5]);dim(Gene_Expression)

Smaple_Summary <- get(load("./Smaples_Summary.Rdata"))
head(Smaple_Summary);dim(Smaple_Summary)


####summary of tissue and cell types
Tissues_s <- Smaple_Summary$Tissue_s
table(Tissues_s)##smaple size per tissue/cell type
length(unique(Tissues_s))##number of unique tissue/cell types


##########################################################################################################
##########################################################################################################
########################To indentify the tissue specific genes############################################


###center and scale all the columens, i.e., normalize across smaples
Gene_Expression_T <- t(Gene_Expression)

head(Gene_Expression_T[,1:5])

###loged and scaled RPKM
Gene_Expression_T_scaled <- log(Gene_Expression_T+0.25)
Gene_Expression_T_scaled <- apply(Gene_Expression_T_scaled,2,scale)

dim(Gene_Expression_T_scaled);class(Gene_Expression_T_scaled);head(Gene_Expression_T_scaled[,1:5])
rownames(Gene_Expression_T_scaled) <- row.names(Gene_Expression_T)
head(Gene_Expression_T_scaled[,1:10])


################################################################################
######################detecting tissue-specific genes using t statistics########
#step 1. group tissues into visual classes, on any subsets
All_samples <- Smaple_Summary$Tissue_s
All_samples_e <- Smaple_Summary$Tissue_e
Samples_sets <- unstack(All_samples,All_samples~as.factor(All_samples_e))
str(Samples_sets)

#step 2. compute the t-statisitcs while correcting for study, age ans sex covariables 
Category <- names(Samples_sets)
tissues_unique <- lapply(Samples_sets,unique)
str(tissues_unique)

##before run, you should change the output dir

for(i in 1:length(Samples_sets)){
  
  Tissue_class <- Category[i]
  No_tissues_category <- length(tissues_unique[[i]])
  ###only one tissue in the category 
  if(No_tissues_category==1){
    
    X <- Gene_Expression_T_scaled
    Sample_class <- ifelse(Smaple_Summary$Tissue_visa%in%Tissue_class,1,-1)
    Sex <- ifelse(Smaple_Summary$Sex%in%c("M"),1,-1)
    Study <- Smaple_Summary$Study
    Age <- Smaple_Summary$Age
    
    Gene_number <- dim(X)[1]
    Myres <- array(NA,dim = c(Gene_number,2))
    
    ##for loop to compute the t-statistics for each  gene in a tissue
    for(i in 1:Gene_number){
      Y <- as.numeric(X[i,])
      myfit <- lm(Y~Sample_class+factor(Study)+factor(Sex)+factor(Age))
      a <- summary(myfit)
      t <- coef(a)["Sample_class",c("t value","Pr(>|t|)")]
      Myres[i,c(1,2)] <- t
      
    }
    row.names(Myres) <- row.names(X)
    write.table(Myres,file = paste("Output_path",Tissue_class,"_t.txt",sep = ""),append = F,quote = F,row.names = T,col.names = F)
  }
  
  ##when tissues belonging to "Other" group
  if(Tissue_class=="Other"){
    for(k in 1:length(tissues_unique[["Other"]])){
      Other_tissues <- tissues_unique[["Other"]]
      X <- Gene_Expression_T_scaled
      Sample_class <- ifelse(Smaple_Summary$Tissue_s%in%Other_tissues[k],1,-1)
      Sex <- ifelse(Smaple_Summary$Sex%in%c("M"),1,-1)
      Study <- Smaple_Summary$Study
      Age <- Smaple_Summary$Age
      
      Gene_number <- dim(X)[1]
      Myres <- array(NA,dim = c(Gene_number,2))
      
      ##for loop to compute the t-statistics for each  gene in rumen
      for(i in 1:Gene_number){
        Y <- as.numeric(X[i,])
        myfit <- lm(Y~Sample_class+factor(Study)+factor(Sex)+factor(Age))
        a <- summary(myfit)
        t <- coef(a)["Sample_class",c("t value","Pr(>|t|)")]
        Myres[i,c(1,2)] <- t
      }
      row.names(Myres) <- row.names(X)
      write.table(Myres,file = paste("Output_path",Other_tissues[k],"_t.txt",sep = ""),append = F,quote = F,row.names = T,col.names = F)
    }
  }
  
  ##when multiple tissues in the category
  if(Tissue_class!="Other"&No_tissues_category>1){
    Tissue_analyzed <- tissues_unique[[i]]
    
    for(j in 1:length(Tissue_analyzed)){
      target_tissue <- Tissue_analyzed[j]
      Similar_Tissues <- Tissue_analyzed[!Tissue_analyzed%in%target_tissue]
      
      X <- Gene_Expression_T_scaled
      
      index <- which(Smaple_Summary$Tissue_s%in%Similar_Tissues)
      
      X <- X[,-index]
      Smaple_Summary_mammary <- Smaple_Summary[-index,]
      
      Sample_class <- ifelse(Smaple_Summary_mammary$Tissue_s%in%target_tissue,1,-1)
      Sex <- ifelse(Smaple_Summary_mammary$Sex%in%c("M"),1,-1)
      Study <- Smaple_Summary_mammary$Study
      Age <- Smaple_Summary_mammary$Age
      
      Gene_number <- dim(X)[1]
      Myres <- array(NA,dim = c(Gene_number,2))
      
      ##for loop to compute the t-statistics for each  gene in a tissue
      for(i in 1:Gene_number){
        Y <- as.numeric(X[i,])
        myfit <- lm(Y~Sample_class+factor(Study)+factor(Sex)+factor(Age))
        a <- summary(myfit)
        t <- coef(a)["Sample_class",c("t value","Pr(>|t|)")]
        Myres[i,c(1,2)] <- t
      }
      
      head(Myres)
      row.names(Myres) <- row.names(RNA_seq)
      write.table(Myres,file = paste("Output_path",target_tissue,"_t.txt",sep = ""),append = F,quote = F,row.names = T,col.names = F)
    }
  }
}


