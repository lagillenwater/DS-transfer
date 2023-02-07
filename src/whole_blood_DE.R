library(tidyverse)
library(ggplot2)
library(limma)


comorbidities <- read.delim("./data/P4C_Comorbidity_020921.tsv", skip = 1)
head(comorbidities)

metadata <-read.delim("./data/P4C_metadata_021921_Costello.txt")
head(metadata)

comorbidities <- comorbidities %>%
    select(RecordID, Condition, HasCondition) %>%
    pivot_wider(names_from = Condition, values_from = HasCondition)




load("./data/DS_transcription_profiles_wide.RData")

ds_wide <- ds_wide %>%
    distinct(LabID, .keep_all = T) # remove duplicate names


## variance filter function
varFilter <- function(gene_expression,threshold) {
    #variance filtere
    variance <- apply(gene_expression[,2:ncol(gene_expression)],2,var)
                                        # numebr of genes in DS cohort with VAR > .1
    filtered <- gene_expression %>%
        dplyr::select(LabID,names(variance[(variance > threshold)]))
    return(filtered)
}



## organize meta
meta <- metadata %>%    inner_join(comorbidities, by = "RecordID")
ds <- varFilter(ds_wide, .1)

meta <- meta %>% filter(LabID %in% ds$LabID) %>% column_to_rownames('LabID') 
meta$numeric_karyotype <- ifelse(meta$Karyotype == "T21", 1, 0)

ds <- ds %>% filter(LabID %in% rownames(meta))  %>% arrange(match(LabID, rownames(meta)))%>% column_to_rownames('LabID') %>% t() %>% as.data.frame()
ds <- log(ds + 1e-10)
##ds <- cbind(ds, karyotype = meta$numeric_karyotype)

design <- model.matrix(~0 + Karyotype +Sex+Age_at_visit, meta)
colnames(design)[1:2] <- c('D21', 'T21')
contrasts <- makeContrasts(T21vsD21 = T21-D21, levels = design)
fit <- eBayes(contrasts.fit(lmFit(ds, design), contrasts))
results <- topTable(fit, adjust = 'fdr', number = nrow(ds) + 1)
results <- results[!is.na(results$logFC), ]
write.csv(results, "whole_blood_DE_results.csv")


T21_meta <- meta %>% filter(Karyotype == "T21")
T21_meta$Depression[is.na(T21_meta$Depression)] <- 0

design <- model.matrix(~0 + Depression +Sex+Age_at_visit, T21_meta)
T21_ds <- ds %>% select(rownames(design))
contrasts <- makeContrasts(Depression, levels = design)
fit <- eBayes(contrasts.fit(lmFit(T21_ds, design), contrasts))
results <- topTable(fit, adjust = 'fdr', number = nrow(ds) + 1)
results <- results[!is.na(results$logFC), ]
write.csv(results, "whole_blood_DE_results.csv")
