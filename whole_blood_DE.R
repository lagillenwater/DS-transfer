library(tidyverse)
library(ggplot2)
library(DESeq2)


comorbidities <- read.delim("./data/P4C_Comorbidity_020921.tsv", skip = 1)
head(comorbidities)

metadata <-read.delim("./data/P4C_metadata_021921_Costello.txt")
head(metadata)

comorbidities <- comorbidities %>%
    select(RecordID, Condition, HasCondition) %>%
    pivot_wider(names_from = Condition, values_from = HasCondition)


meta <- metadata %>%
    inner_join(comorbidities, by = "RecordID")

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

ds <- varFilter(ds_wide, .1)
ds <- ds %>%
    column_to_rownames('LabID')

ds <- as.data.frame(t(ds))

# organize meta

ds <- DESeqDataSetFromMatrix(countData = ds,colData = ,design= ~ batch + condition)

