#### This program just reshapes the gene expression counts from long to wide

library(tidyverse)

## load expression data
expression <- read.delim("./data/HTP_WholeBlood_RNAseq_Counts_Synapse.txt")


## pivot wider. Have to include the unique identifier for row numbers to avoid any errors. 
htp_expr <- expression %>%
    select(LabID, Gene_name,EnsemblID,Chr,Value) %>%
    group_by(Gene_name) %>%
    mutate(row = row_number())%>%
    pivot_wider(names_from = LabID, id_cols = c(Gene_name, Chr,EnsemblID), values_from = Value) 

save(htp_expr, file = "./data/HTP_transcription_counts_wide.Rdata")

