## Sample data processing

library(tidyverse)

## Create smaller sample datasets for testing

## Prostate Cancer Data
meta <- read.csv("./data/tcga/PRAD_metadata.csv")
expr <- read.csv("./data/tcga/PRAD_expression.csv")

## Identify  ID's with primary tumor and normal tissue samples
tcga_ids <- meta %>% select(tcga.cgc_case_id,tcga.cgc_sample_sample_type)  %>% filter(tcga.cgc_sample_sample_type == "Solid Tissue Normal") %>% .$tcga.cgc_case_id

## sample ids are different from tcga ids
sample_ids <- meta %>% filter(tcga.cgc_case_id %in% tcga_ids) %>% .$X

## filter expression  data and metadata by those ids
sample_expression <- expr %>% filter(X %in% sample_ids)
sample_meta <- meta %>% filter(tcga.cgc_case_id %in% tcga_ids)

    
## Apply variance filter to genes
gene_vars <- sample_expression %>%
    summarise(across(where(is.numeric), var))

high_variance_genes <- gene_vars %>%
    pivot_longer(cols = everything(), names_to = "gene",values_to = "variance") %>%
    arrange(desc(variance)) %>%
    filter(variance > 1)

sample_expression_filtered <- sample_expression %>%
    dplyr::select(c('X', high_variance_genes$gene))

## log transform
sample_expression_filtered[,2:ncol(sample_expression_filtered)] <- log2(sample_expression_filtered[,2:ncol(sample_expression_filtered)]+1)

## Convert Ensembl IDs to Gene Names
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- names(sample_expression_filtered)[2:ncol(sample_expression_filtered)]
alternate <- gsub(".*[.]","", genes)
genes <- gsub("[.].*","", genes)
G_list <- getBM(filters= "ensembl_gene_id",attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
G_list <- G_list %>% filter(hgnc_symbol != "")

names(sample_expression_filtered) <- gsub("[.].*","", names(sample_expression_filtered))
sample_expression_filtered <- sample_expression_filtered %>%
    dplyr::select(c("X", G_list$ensembl_gene_id)) %>%
    subset(.,select = which(!duplicated(names(.))))
names(sample_expression_filtered)[2:ncol(sample_expression_filtered)] <- G_list$hgnc_symbol


## DEG
library(limma)
sample_expression_filtered <- sample_expression_filtered %>%
    column_to_rownames('X') %>%
    t()
design <- model.matrix(~0 + as.factor(tcga.cgc_sample_sample_type), data = sample_meta)
colnames(design) <- c("Primary_Tumor", "Solid_Tissue_Normal")
contrast <- makeContrasts(constrast = Primary_Tumor -  Solid_Tissue_Normal, levels = design)
fit <- eBayes(contrasts.fit(lmFit(sample_expression_filtered, design), contrast))
deg_res <- topTable(fit, adjust = 'fdr', number = 13184, p.value = .05)
deg_res <- deg_res %>% arrange(desc(logFC))
write.csv(deg_res,"./results/tcga/deg/prostate_deg.csv", row.names = F)



## PCA exploration
library(FactoMineR)

data <- sample_expression_filtered %>%
    inner_join(sample_meta %>% dplyr::select(X,tcga.cgc_sample_sample_type, tcga.cgc_case_id, tcga.cgc_case_age_at_diagnosis ), by = 'X') %>%
    filter(!(tcga.cgc_case_id %in% c("TCGA-HC-8258", "TCGA-HC-7740"))) %>%
    column_to_rownames('X') 

res.pca <- PCA(data, quanti.sup = c("tcga.cgc_case_age_at_diagnosis" ), quali.sup = c("tcga.cgc_sample_sample_type", "tcga.cgc_case_id"), graph = F, ncp = 100)
plot.PCA(res.pca, label = "none")

PC1 <- res.pca$ind$coord[,1]
PC2 <- res.pca$ind$coord[,2]
labs <- data$tcga.cgc_case_id
PCs <- data.frame(cbind(PC1,PC2))

# Just showing the individual samples...
p1 <- ggplot(PCs, aes(PC1,PC2, color = data$tcga.cgc_sample_sample_type)) + 
    geom_point() +
    theme_minimal() +
    theme(legend.position = "bottom")



p2 <- ggplot(PCs, aes(PC1,PC2, color = data$tcga.cgc_case_age_at_diagnosis)) + 
    geom_point() +
    theme_minimal() +
    theme(legend.position = "bottom")


library(gridExtra)

grid.arrange(p1,p2)



#### Old Code

## Considering filtering genes by named pathways. After consulation with Mike about the potential for pathways including genes not relevant to prostate cancer, I decided to test empirically for differentially expressed genes. 

## ## Filter genes by pathways
## library('biomaRt')
## mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
## genes <- names(sample_expression)[2:ncol(sample_expression)]
## genes <- gsub("[.].*","", genes)

## G_list <- getBM(filters= "ensembl_gene_id",
##attributes= c("ensembl_gene_id","hgnc_symbol"),
##values=genes,mart= mart)

## names(sample_expression) <- gsub("[.].*","", names(sample_expression))


## ## Prostate pathways
## pro <- read.delim("../data/Gene Pathways/KEGG_PROSTATE_CANCER.v2023.1.Hs.tsv")
## pro <- pro[17,2]
## pro <- unlist(strsplit(unlist(pro), ","))

## ## Bladder pathways
## bla <- read.delim("../data/Gene Pathways/KEGG_BLADDER_CANCER.v2023.1.Hs.tsv")
## bla <- bla[17,2]
## bla <- unlist(strsplit(unlist(bla), ","))


## ## Housekeeping genes
## h <- read.delim("../data/Gene Pathways/HSIAO_HOUSEKEEPING_GENES.v2023.1.Hs.tsv")
## h <- h[17,2]
## h <- unlist(strsplit(unlist(h), ","))

## filter_genes <- unique(c(pro,bla,h))
## filter_genes <- filter_genes[filter_genes != ""]

## G_list <- G_list %>%
##     filter(hgnc_symbol %in% filter_genes)

## sample_expression_filtered <- sample_expression %>%
##     dplyr::select(c("X", G_list$ensembl_gene_id)) 


## names(sample_expression_filtered)[2:ncol(sample_expression_filtered)] <- G_list$hgnc_symbol

## sample_expression_filtered <- subset(sample_expression_filtered, select = which(!duplicated(names(sample_expression_filtered))))
