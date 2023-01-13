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

## This function finds the total number of individuals in the study, including controls. 
findVariable_controls <- function(metadata_table, variable) {
    if(any(grepl(variable, metadata_table$V1, ignore.case = T))) {
        metadata  <- metadata_table$V1
    } else {
        metadata <- character(0)
    }
           
    return(metadata)
}

## findVariableWrapper is a wrapper function for applying the findVariable over a list of metadata data frames.
findVariableWrapper <- function(metadata, variable, FUN = findVariable() ) {
    variable_rows <- lapply(metadata, function(x) {FUN(x, variable)})
    variable_rows <- variable_rows[lapply(variable_rows,length)>0]
    return(variable_rows)
}


## variableTable is a function for creating a table of counts of metadata by vector of variables
variableTable <- function(metadata, variables, FUN = findVariable()) {
    variable_counts <- lapply(variables, function(x) {sum(unlist(lapply(findVariableWrapper(metadata,x, FUN), length)))})
    variable_table <- data.frame(variable = variables, count = unlist(variable_counts))
    return(variable_table)
}

## find study is a helper function to print information on a particular study
findStudy <- function(srp,metadata_dir = "./recount_metadata/", metadata=metadata) {
    metadata_files <- paste0(metadata_dir, list.files(metadata_dir))
    print(metadata[[grep(srp,metadata_files)]])
}





metadata <- readMetadataWrapper("./recount_metadata/")

tissues <- c("blood", "brain", "breast", "fibroblast", "lymphoblast", "bladder", "colon", "heart", "liver", "muscle","prostate", "skin", "pancreas", "lung", 'testis', 'spleen', 'thyroid', 'ovary', 'esophagus', 'kidney', 'salivary', 'small intestine', 'stomach', 'uterus', 'vagina', 'bone', 'scalp')
tissues <- tissues[order(tissues)]
tissue_table <- variableTable(metadata, variables = tissues, FUN = findVariable_controls)
tissue_table

library(ggplot2)
ggplot(tissue_table, aes(x = variable, y = count)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 60,vjust = 0.5))


conditions <- c("alzh", 'autism', 'epilep', 'leuk', 'autoimmune', 'alopecia', 'arthritis', 'celiac', 'diabet', 'sleep', 'heart', 'thyroidism', 'depression',  'obes','seizure',  'T21')
conditions <- c("leukemia", "lin-", "aml", "mll")
conditions <- c("autoimmune","auto-immune", "immune")
conditions <- c("thyroidism")
conditions <- c("diabetes", "islet")
conditions <- conditions[order(conditions)]
condition_table <- variableTable(metadata,conditions, FUN = findVariable)
condition_table_controls <- variableTable(metadata,conditions, FUN = findVariable_controls)
condition_table
condition_table_controls




ggplot(condition_table, aes(x = variable, y = count)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 60,vjust = 0.5))
