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
    metadata_vector <- metadata_table[grepl
    (variable, metadata_table$V1, ignore.case = TRUE),]
    return(metadata_vector)
}


## excludeVariable is a function for excluding metadata with a full term that differs from the partial term of interest. For example, when looking for 'blood' samples exclude 'blood vessel'
excludeVariable <- function(metadata_vector,variable) {
    metadata_vector <- metadata_vector[!(grepl(variable, metadata_vector, ignore.case = TRUE))]
    return(metadata_vector)
}
                                 
## This function finds the total number of individuals in the study, including controls. 
findVariable_controls <- function(metadata_table, vaqriable) {
    if(any(grepl(variable, metadata_table$V1, ignore.case = TRUE))) {
        metadata_vector  <- metadata_table$V1
    } else {
        metadata_vector <- character(0)
    }
    return(metadata_vector)
}

## containsVariables is a function that takes a vector of variables and identifies those containing similar strings. For example, "blood" and "blood vessel" should identify distinct metadata. The output of this function can feed into the "exclude variable" function.
containsVariables <- function(variable, variables) {
    overlapping_variables <- variables[grepl(variable, variables)]
    overlapping_variables <- overlapping_variables[overlapping_variables != variable]
    return(overlapping_variables)
}

## findVariableWrapper is a wrapper function for applying the findVariable over a list of metadata data frames.
findVariableWrapper <- function(metadata, variable, FUN = findVariable ) {
    variable_rows <- lapply(metadata, function(x) {
        FUN(metadata_table = x, variable)
    })
    variable_rows <- variable_rows[lapply(variable_rows,length)>0]
    return(variable_rows)
}


## excludeVariableWrapper is a wrapper function for applying the findVariable over a list of metadata data frames.
excludeVariableWrapper <- function(metadata_list, variable ) {
    variable_rows <- lapply(metadata_list, function(x) excludeVariable(x,variable))
    variable_rows <- variable_rows[lapply(variable_rows,length)>0]
    variable_rows <- unlist(variable_rows)
    return(variable_rows)
}

## variableTable is a function for creating a table of counts of metadata by vector of variables
variableTable <- function(metadata, variables, FUN = findVariable) {
    variable_counts <- lapply(variables, function(x) {
        res <- findVariableWrapper(metadata,x, FUN)
        similar_variables <- containsVariables(x,variables)
        res <- unlist(res)
        if(length(similar_variables) >0) {
           for( y in similar_variables){
                tmp <- unlist(excludeVariableWrapper(res,y))
                tmp <- unlist(tmp)
                res <- res[!(res %in% intersect(res,tmp))]
            }
        }
        return(length(res))
    })
    variable_table <- data.frame(variable = variables, count = as.numeric(unlist(variable_counts)))
    return(variable_table)
}

## find study is a helper function to print information on a particular study

findStudy <- function(srp,metadata_dir = "./recount_metadata/", metadata=metadata) {
    metadata_files <- paste0(metadata_dir, list.files(metadata_dir))
    return(metadata[[grep(srp,metadata_files)]])
}

## function for counting all the metadata
countMetadata <- function(metadata) {Figure 1: Barplot depicting the count of tissue type within the samples that had the T21 kary- otype based on the key term search.
    count <- sum(unlist(lapply(metadata,nrow)))
    return(count)
}

## variableBarPlots is a function for creating barplots based on the output of variableTable
library(ggplot2)
variableBarPlots <- function(variable_table) {
    p1 <- ggplot(variable_table, aes(x = variable, y = count)) +
        geom_bar(stat = "identity") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 12),
             axis.text.y = element_text(size = 12),
             legend.position = "none") +
        xlab("")
    return(p1)
}


### The analysis portion

metadata <- readMetadataWrapper("./recount_metadata/")


count <- countMetadata(metadata)

tissues <- c( "brain", "breast", "fibroblast", "lymphoblast", "bladder", "colon", "heart", "liver", "muscle","prostate", "skin", "pancreas", "lung", 'testis', 'spleen', 'thyroid', 'ovary', 'esophagus', 'kidney', 'salivary', 'small intestine', 'stomach', 'uterus', 'vagina',  'scalp', 'pbmc', 'plasma', 'adipose','adrenal','blood vessel', 'bone marrow', 'brain reference', 'cervix','colon','colorectal','epithelium','esophagus','glioblastoma','head and neck', 'large intestine','melanoma','muscle','nerve','ovary','pituitary','soft tissue', 'spleen','thymus','testes','tonsil','umbilical cord',  "whole blood", "PBMC","fibroblast", "lymphoblast", "blood vessel", "ipsc", "peripheral blood mononuclear cell", "leukocyte", "monocyte", "lymphocyte", "monocyte")


tissues <- unique(tolower(tissues[order(tissues)]))
tissue_table <- variableTable(metadata, variables = tissues, FUN = findVariable)
tissue_table

tissue_table <- tissue_table[tissue_table$variable %in% tissues,]
tissue_table[33, "count"] = tissue_table[33, "count"] + tissue_table[32, "count"]
tissue_table <- tissue_table[-32,]
write.csv(tissue_table, file = "tissue_table.csv", row.names = F, quote = F)
#tissue_table <- rbind(tissue_table, c("not annotated", count - sum(tissue_table$count)))
tissue_plot <- variableBarPlots(tissue_table)

jpeg( file = "tissue_plot.jpeg")
print(tissue_plot)
dev.off()


### Blood breakdwon

## A problem with using these terms is that there are often synonyms to these terms that prevent exact matching. Or include an overlap between the query term and the synonym. A better approach may be to use these synonyms in the search.
## A resource could be the node identifiers from pheknowlator. This includes data like synonyms. 
kg <- read.csv("filtered_identifier.csv")
pbmc <- kg[grepl("peripher", kg$Label, ignore.case = TRUE) | grepl("periph", kg$synonym, ignore.case = TRUE), ]

## Upon searching for terms related to PBMC's, I was able to locate some terms with the fuzzy search. For example, the term "periph" returned 298 hits. Hit #151 was for peripheral blood mononuclear cells. However, term didn't have any synonyms. At least for this example, this wasn't helpful. Maybe worth returning to later with another resource. 

blood <- c( "whole blood", "lymphoblast", "blood vessel", "ipsc", "peripheral blood mononuclear cell", "leukocyte", "monocyte", "lymphocyte", "monocyte", "plasma")
blood <- unique(blood[order(blood)])
blood_meta <- lapply(blood, function(x) findVariableWrapper(metadata, x))
blood_meta <- unlist(blood_meta)
blood_meta <- lapply(blood_meta, as.data.frame) # findVariable converts the data to a character vector while variableTable is looking for a data frame
blood_meta <- lapply(blood_meta,setNames, "V1")
blood_count <- countMetadata(blood_meta)
blood_table <- variableTable(blood_meta, variables = blood, FUN = findVariable)

blood_table <- tissue_table[tissue_table$variable %in% blood,]

blood_plot <- variableBarPlots(blood_table)
jpeg( file = "blood_plot.jpeg")
print(blood_plot)
dev.off()
write.csv(blood_table, file = "blood_table", row.names = F, quote = F)


conditions <- c("alzh", 'autism', 'epilep', 'leuk', 'autoimmune', 'alopecia', 'arthritis', 'celiac', 'diabet', 'sleep', 'heart', 'thyroidism', 'depression',  'obes','seizure',  'T21', 'Trisomy21', 'Down syndrome', "Down's syndrome", "patient21", "subject21","s21","TET21", "chromosome 21", "DS", "trisomic", "trisomy 21")
conditions <- conditions[order(conditions)]
condition_table <- variableTable(metadata,conditions, FUN = findVariable)




whole_meta <- findVariableWrapper(metadata, "whole blood")
whole_meta <- lapply(whole_meta, as.data.frame)
whole_meta <- lapply(whole_meta, setNames, "V1")


condition_table_whole <- variableTable(whole_meta, variables = conditions, FUN=findVariable)






                                        #condition_table_controls <- variableTable(metadata,conditions, FUN = findVariable_controls)


condition_table
condition_table[2, "variable"] <- "Alzheimer's"
condition_table[condition_table$variable == "trisomy 21", "count"] = 133
condition_table[c(9,13,14,15,16), "variable"] <- c("Diabetes", "Epilepsy", "Heart Defects",  "Leukemia", "Obesity")

condition_table <- condition_table[-c(7,10,11,12,17,18,21,22,23,25,27),]

condition_table

condition_plot <- variableBarPlots(condition_table)
write.csv(condition_table, file = "condition_table", row.names = F, quote = F)

jpeg( file = "condition_plot.jpeg")
print(condition_plot)
dev.off()

T21_variables <- c('T21', 'Trisomy21', 'Down syndrome', "Down's syndrome",  "chromosome 21",  "trisomic", "trisomy 21")
T21_meta <- lapply(T21_variables, function(x) findVariableWrapper(metadata, x))
T21_meta <- unlist(T21_meta)

T21_meta <- excludeVariableWrapper(T21_meta, "patient21")
T21_meta <- excludeVariableWrapper(T21_meta, "subject21")
T21_meta <- excludeVariableWrapper(T21_meta, "TET21")
T21_meta <- excludeVariableWrapper(T21_meta, "pet21")
T21_meta <- excludeVariableWrapper(T21_meta, "pat21")
T21_meta <- excludeVariableWrapper(T21_meta, "07T21")
T21_meta <- excludeVariableWrapper(T21_meta, "20T21")
T21_meta <- excludeVariableWrapper(T21_meta, "17T21")
T21_meta <- excludeVariableWrapper(T21_meta, "27T21")
T21_meta <- excludeVariableWrapper(T21_meta, "11T21")
T21_meta <- excludeVariableWrapper(T21_meta, "05T21")
T21_meta <- excludeVariableWrapper(T21_meta, "03T21")
T21_meta <- excludeVariableWrapper(T21_meta, "10T21")
T21_meta <- excludeVariableWrapper(T21_meta, "16T21")
T21_meta <- excludeVariableWrapper(T21_meta, "08T21")
T21_meta <- excludeVariableWrapper(T21_meta, "T21B")
T21_meta <- excludeVariableWrapper(T21_meta, "CT21")
T21_meta <- excludeVariableWrapper(T21_meta, "02T21")
T21_meta <- excludeVariableWrapper(T21_meta, "Alias;;T21")
T21_meta <- unlist(T21_meta)
T21_meta <- lapply(T21_meta, as.data.frame) # findVariable converts the data to a character vector while variableTable is looking for a data frame
T21_meta <- lapply(T21_meta, setNames, "V1")
T21_count <- countMetadata(T21_meta)
tissues <- c(tissues, "T cell")

T21_table <- variableTable(T21_meta, variables = tissues, FUN = findVariable)
T21_table$count <- as.numeric(T21_table$count)

T21_table <- rbind(T21_table, c("not annotated", T21_count - sum(T21_table$count)))
T21_table <- T21_table[T21_table$count > 0,]

T21_plot <- variableBarPlots(T21_table)
jpeg( file = "T21_plot.jpeg")
print(T21_plot)
dev.off()
write.csv(T21_table, file = "T21_table", row.names = F, quote = F)


ggplot(condition_table, aes(x = variable, y = count)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 60,vjust = 0.5))





findVariableWrapper(metadata, "depress")
