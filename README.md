# Cattle-GeneAtlas
Gene expression (FPKM) across 91 tissues and cell types from 723 RNA-seq datasets in cattle. We computed the t-statistic to measure the tissue-specific expression pattern for each gene in each tissue and cell type. Our cattle-GeneAtlas V1 provides a primary source for GWAS interpretation, functional validation, studies of adaptive evolution, domestication and genomic improvement in livestock. For intance, you could conducted a tissue-enrichment analysis (like GO or KEGG enrichment analysis) using Cattle-GeneAtlas to detect the potential tissues or cell types on which  genes of your interest (e.g., from a certain GWAS or evolution study) would act. You also could differentially weight SNPs in genes that are highly specificially expressed in the trait-relevant tissues in the prediction models (like SSGP, multi-BLUP or BayesRC) to improve genomic prediction for that particular trait, especially in the multi-breeds or over-generation situations.

# Cattle_GeneAtlas_Expression_Summary_Clustering.R
This R code can be used to explore and plot the number of expressed genes and clustering samples using both PCA and t-SNE approaches

