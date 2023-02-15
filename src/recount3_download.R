## This is an Rscript to be run as a command line program for downloading and processing specific gene expression files 

## Checking for libraries.
## TODO rewrite to work run in packages format. 

if (!requireNamespace("BiocManager", quietly = TRUE)) { 
    install.packages("BiocManager")
}

if (!requireNamespace("optparse", quietly = TRUE)) { 
    install.packages("optparse")
}
if (!requireNamespace("recount3", quietly = TRUE)) { 
    BiocManager::install("LieberInstitute/recount3") # install recount3
    BiocManager::install("recount")

}
if (!requireNamespace("recount", quietly = TRUE)) { 
    BiocManager::install("LieberInstitute/recount3") # install recount3
    BiocManager::install("recount")

}
if (!requireNamespace("tidyverse", quietly = TRUE)) { 
    install.packages("tidyverse")
}

## Load required libraries 

library("optparse")
library(recount3) 
library(recount)
library(tidyverse)

option_list = list(
    make_option(c("-f", "--file_source"), type="character", default=NULL, 
              help="transcriptomic database ('sra', 'gtex', 'tcga')", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt", 
                help="output file name [default= %default]", metavar="character"),
    make_option(c("-p", "--project"), type="character", default=NULL, 
                help="project number or title", metavar="character"),
    make_option(c("-s", "--species"), type="character", default='human', 
                help="'human' or 'mouse' [default='human']", metavar="character")    
   
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## TODO ### Check for NULLs

human_projects <- available_projects(organism = opt$species) # load all human projects

human_source <- human_projects %>% filter(file_source == opt$file_source)

## TODO Add fuzzy match
project <- human_source %>% filter(project == opt$project)

print(project)

## TOOD Add writing of file
