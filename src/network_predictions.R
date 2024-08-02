
####This  program uses gene regulatory networks to predict gene expression based on transcription factor profiles.

### libraries
library(tidyverse)
library(caret)
library(ggplot2)
library(foreach)
library(doParallel)



load("./data/DS_transcription_profiles_wide.RData") # gene expression profiles in a wide format, samples x genes. Data comes fromo whole blood samples. 
ds_wide <- ds_wide %>%
    distinct(LabID, .keep_all = TRUE) %>%
    as.data.frame        # For some reason the samples are duplicated. Need to preprocess before data can be command line input.
names(ds_wide) <- make.names(names(ds_wide))
metadata <-read.delim("./data/P4C_metadata_021921_Costello.txt") # metadata for gene expression data. Include Sex, Age, Karyotype, Event_name, Age, BMI, and samle source
load("../graphTransform/data/GTEx_PANDA_tissues.RData") # tissue specific gene regulatory networks created by Sonawane et al.



## variance filter function
varFilter <- function(gene_expression,threshold) {
    #variance filtere
    variance <- apply(gene_expression[,2:ncol(gene_expression)],2,var)
                                        # numebr of genes in DS cohort with VAR > .1
    filtered <- gene_expression %>%
        dplyr::select(names(variance[(variance > threshold)]))
    return(filtered)
}



### Function to convert ensembl IDs in the network into gene names
convertENSG <- function(network,gene_meta) {
    gene_meta$Symbol <- make.names(gene_meta$Symbol)
    network <- network %>% inner_join(gene_meta, by = c("Gene" = "Name"))
    names(network)[c(1,4)] <- c("source", "target")
    return(network)
}



# Filter the network to only include tuples that coNtai genes from the gene list in the source and target
filterNetwork <- function(network, gene_vector) {
    network <- network %>%
        filter(source %in% gene_vector & target %in% gene_vector)
    return(network)
}    


# Add the tissue specific weight to the network. 
tissueNetwork <- function(network,tissue_weights, tissue) {
    network <- cbind(network, weight = tissue_weights[, tissue])
    return(network)
}

### find the TFs for a gene
findTF <- function(network, target_gene) {
    TF <- network %>% filter(target == target_gene) %>% .$source
    return(TF)
}


## predict gene expression values using incoming weights and expression values.
predictExpression <- function(network,  target_gene, gene_expression) {
    target_network <- network %>% filter(target == target_gene) %>% filter(weight > 0)
    if(nrow(target_network) <= 1) {
        prediction <- rep(NA, nrow(gene_expression))
        } else {
        
   # target_network <- target_network[1:10,]
            TF_expression_profiles <- gene_expression %>%
                dplyr::select(target_gene, target_network$source)

                                        # define training control
            train_control <- trainControl(method="cv", number=5)
            model_formula <- as.formula(paste(target_gene, '~.'))
            model <- train(model_formula , data = TF_expression_profiles, trControl=train_control, method="lm")
            prediction <- predict(model, TF_expression_profiles)
            }
    return(prediction)
}

## compare correlation of predictions between karyotype
predictionAccuracy <- function(network, target_gene, gene_expression_1, gene_expression_2) {
    
    prediction_1 <- suppressWarnings(predictExpression(network, target_gene, gene_expression_1))
    prediction_2  <- suppressWarnings(predictExpression(network, target_gene, gene_expression_2))
    cor_1 <- cor(prediction_1, gene_expression_1[, target_gene])
    cor_2 <- cor(prediction_2, gene_expression_2[, target_gene])
    return(c(cor_1,cor_2))
} 



### wrapper to generate predicitons over all genes
predictionWrapper <- function(network,tissue_weights,gene_meta, gene_expression,expression_metadata, variance_threshold = .1, tissue = "Whole_blood", chr21) {
    network <- convertENSG(network, gene_meta)
    network <- tissueNetwork(network, tissue_weights,tissue)
    
    gene_expression_1 <- gene_expression %>%
        filter(LabID %in% (expression_metadata %>% filter(Karyotype == "T21") %>% .$LabID))
    gene_expression_2 <- gene_expression %>%
        filter(LabID %in% (expression_metadata %>% filter(Karyotype == "Control") %>% .$LabID))
    
    gene_expression_1 <- varFilter(gene_expression_1, threshold = variance_threshold)
    gene_expression_2 <- varFilter(gene_expression_2, threshold = variance_threshold)

    genes_to_keep <- intersect(intersect(names(gene_expression_1), names(gene_expression_2)), network$target)
    genes_to_keep <- c(genes_to_keep[genes_to_keep %in% chr21$symbol], sample(genes_to_keep,(1000-108)))
    genes_to_keep <- unique(genes_to_keep)
    

    gene_expression_1 <- gene_expression_1 %>% dplyr::select(genes_to_keep)
    gene_expression_2 <- gene_expression_2 %>% dplyr::select(genes_to_keep)

    network <- filterNetwork(network, as.character(genes_to_keep))

    targets <- as.character(unique(network$target))
    targets <- targets[targets %in% genes_to_keep]
    length_targets <- length(targets)
    predicted_accuracy <- list()

    #numCores <- detectCores()

    #registerDoParallel(numCores-2)
    for(x in  1:length_targets) {
        print(x/length_targets)
        predicted_accuracy[[x]] <- predictionAccuracy(network, targets[x], gene_expression_1, gene_expression_2)
    }
    predicted_accuracy <- data.frame(matrix(unlist(predicted_accuracy), nrow=length(targets), byrow=T))
    rownames(predicted_accuracy) <- targets
    names(predicted_accuracy) <- c("T21", "Control")
    return(predicted_accuracy)
}

library(openxlsx)
chr21 <- read.delim("./data/non_alt_loci_set_chr_21.txt")
chr21 <- chr21 %>% filter(e.symbol %in% names(ds_T21))



prediction_accuracy <- predictionWrapper(network = edges,
                                          tissue_weights = net,
                                          gene_meta = genes,
                                           gene_expression = ds_wide,
                                           expression_metadata = metadata,
                                           variance_threshold = .1,
                                         tissue = "Whole_blood",
                                         chr21 = chr21)

toplot <- prediction_accuracy %>%
    rownames_to_column() %>%
    pivot_longer(cols = c("T21", "Control"),names_to = "Karyotype", values_to = "prediction")

toplot$chr21_gene <- ifelse(toplot$rowname %in% chr21$symbol, "yes", "no")
prediction_accuracy$chr21_gene <- ifelse(rownames(prediction_accuracy) %in% chr21$symbol, "yes", "no")

ggplot(prediction_accuracy, aes(x =T21, y = Control, color = chr21_gene)) +
    geom_point(size = 2) +
    geom_abline(intercept = 0, linetype='dashed', col = 'red', size = 1.5)+
    theme_classic() +
    xlab("T21 Gene Expression Prediction Accuracy") +
    ylab("D21 Gene Expression Prediction Accuracy") +
    scale_color_manual(values = c("yes" = "black", "no" = "gray60")) +
    theme(legend.position = "none") +
    xlim(c(0,1)) + ylim(c(0,1))


T21_prediction <- prediction_accuracy %>% filter(chr21_gene == "yes")



rand <- rnorm(nrow(prediction_accuracy), mean = .1, sd = .02)
prediction_accuracy$rand <- ifelse((prediction_accuracy$T21+ rand) >= 1, 1ggplot(toplot, aes(x = factor(Karyotype),y = prediction, color = factor(Karyotype)))+
    geom_boxplot(outlier.alpha = 0)+
    geom_jitter() +
    theme_classic() +
    xlab("") +
    theme(legend.position = "none") +
    ylab("Gene Expression Prediction Accuracy")

t.test(prediction ~ Karyotype, data =toplot)
