### This is a script for reading in the metadata preprocessed in the recount3_tissue.R. That script took all the metadata files and extracted the main variable of interest: "sra.sample_attributes". The txt files containing that information is located in the directory recount3_metadata.



## readMetadata is a script for reading in the text metadata file. It takes as an input the metadata file.
readMetadata <- function(metadata_file) {
    metadata <- read.delim(metadata_file, header = F)
    return(metadata)
}


## readMetadataWrapper is a function for applying the metadata wrapper over all the files in the recount3_metadata directory.
readMetadataWrapper <- function(metadata_dir) {
    metadata_files <- paste0(metadata_dir, list.files(metadata_dir))
    metadata <- lapply(metadata_files, function(x) {
        tryCatch({readMetadata(x)}, error = function(e) {})})      
    return(metadata)
}


##### NB: could improe readMetadataWrapper or readMetadata by converting the text file to a vectore instead of a data frame. Would save space

## 7510 studies have metadata

## findVariable is a function for finding the rows that contain the keyword of interest in the metadata tables.
## Takes as an input the metadata table to search and the variable of interest. 
findVariable <- function(metadata_table, variable) {
    metadata <- metadata_table$V1[grepl(variable, metadata_table$V1, ignore.case = T)]
    return(metadata)
}

## findVariableWrapper is a wrapper function for applying the findVariable over a list of metadata data frames.
findVariableWrapper <- function(metadata, variable) {
    variable_rows <- lapply(metadata, function(x) {findVariable(x, variable)})
    variable_rows <- variable_rows[lapply(variable_rows,length)>0]
    return(variable_rows)
}


## variableTable is a function for creating a table of counts of metadata by vector of variables
variableTable <- function(metadata, variables) {
    variable_counts <- lapply(variables, function(x) {sum(unlist(lapply(findVariableWrapper(metadata,x), length)))})
    variable_table <- data.frame(variable = variables, count = unlist(variable_counts))
    return(variable_table)
}


tissues <- c("blood", "brain", "breast", "fibroblast", "lymphoblast", "bladder", "colon", "heart", "liver", "muscle","prostate", "skin", "pancreas", "lung", 'testis', 'spleen', 'thyroid', 'ovary', 'esophagus', 'kidney', 'salivary', 'small intestine', 'stomach', 'uterus', 'vagina', 'bone', 'scalp')
tissues <- tissues[order(tissues)]
tissue_table <- variableTable(metadata, tissues)

library(ggplot2)
ggplot(tissue_table, aes(x = variable, y = count)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 60,vjust = 0.5))


conditions <- c("alzheimer", 'autism', 'epilepsy', 'leukemia', 'immune', 'alopecia', 'arthritis', 'celiac', 'diabet', 'sleep', 'heart', 'thyroidism', 'depression',  'obes','seizure',  'T21')
conditions <- conditions[order(conditions)]
condition_table <- variableTable(metadata,conditions)


ggplot(condition_table, aes(x = variable, y = count)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 60,vjust = 0.5))
