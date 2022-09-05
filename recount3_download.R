                                        # This is a program for downloading the preliminary data from recount3 into the directory

if (!requireNamespace("BiocManager", quietly = TRUE)) { 
    install.packages("BiocManager")
}

## Check that you have a valid Bioconductor installation
BiocManager::valid()

BiocManager::install("LieberInstitute/recount3") # install recount3
BiocManager::install("recount")

# Load required libraries
library(recount3) # load the recount3 library
library(recount)
library(tidyverse)

human_projects <- available_projects() # load all human projects

dim(human_projects)

head(human_projects)

T21_studies <- c("SRP039348","SRP188973", "SRP188969",  "SRP186520", "SRP017123") # T21 projects




# create list for storing counts
T21_data = list()
T21_metadata = list()

for(i  in 1:length(T21_studies)){
    print(i)
  # get project info for a specific project
  proj_info <- subset(
    human_projects,
    project == T21_studies[i] #& project_type == "data_sources"
    )
  # select 1 project
  rse_gene = create_rse(proj_info)
  assay(rse_gene, "counts") = transform_counts(rse_gene)
  assays(rse_gene)$RPKM = recount::getRPKM(rse_gene)
  metadata =  tryCatch({  
  # metadata 
      expand_sra_attributes(rse_gene)
  }, error = function(e) {})
  if(is.null(metadata)) {
        print(paste("skipping", i))
        next
  }
  metadata = colData(metadata)
  metadata = metadata[,grep("sra_attribute", names(metadata))] 
  # store expression data
  T21_data[[i]] =as.data.frame(t(assays(rse_gene)$RPKM))
  # store the metadata
  T21_metadata[[i]] = metadata
}

# combine list object
T21_expression = do.call(rbind,T21_data)


T21_metadata[[1]]$sra_attribute.karyotype <- ifelse(T21_metadata[[1]]$`sra_attribute.disease/status` == "Trisomy 21", "T21", "D21")
T21_metadata[[4]]$sra_attribute.karyotype <- ifelse(T21_metadata[[4]]$sra_attribute.condition == "Trisomic", "T21", "D21")

T21_meta = lapply(T21_metadata, function(x) x[, c( "sra_attribute.karyotype", "sra_attribute.source_name")])
T21_meta = do.call(rbind, T21_meta)


T21_individuals <- rownames(T21_meta[T21_meta$sra_attribute.karyotype == 'T21',]) # 23 people with DS karyotype

T21_ge <- T21_expression[T21_individuals,]


### Potentially follow up here with GTEX samples, the "healthy" ones. 
gtex <- human_projects %>%
    filter(file_source == 'gtex')

blood <- gtex %>% # just use blood for this first analysis
    filter(project == "BLOOD")

# create list for storing counts
gtex_data = list()
gtex_metadata = list()

for(i  in blood$project){
    print(i)
  # get project info for a specific project
  proj_info <- subset(
    human_projects,
    project == i #& project_type == "data_sources"
    )
  # select 1 project
  rse_gene = create_rse(proj_info)
  assay(rse_gene, "counts") = transform_counts(rse_gene)
  assays(rse_gene)$RPKM = recount::getRPKM(rse_gene)
  metadata =  tryCatch({  
  # metadata 
      expand_sra_attributes(rse_gene)
  }, error = function(e) {})
  if(is.null(metadata)) {
        print(paste("skipping", i))
        next
  }
  metadata = colData(metadata)
  #metadata = metadata[,grep("sra_attribute", names(metadata))] 
  # store expression data
  gtex_data[[i]] =as.data.frame(t(assays(rse_gene)$RPKM))
  # store the metadata
  gtex_metadata[[i]] = metadata
}

gtex_expression <- as.data.frame(gtex_data)
gtex_meta <- as.data.frame(gtex_metadata)


#restricn data to JAK/STAT pathway
library(msigdbr)
all_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
all_gene_sets <- as.data.frame(all_gene_sets)
jak_stat <- all_gene_sets[all_gene_sets$gs_name == "HALLMARK_IL6_JAK_STAT3_SIGNALING",]


T21_jak <- T21_ge[, grepl(paste(jak_stat$ensembl_gene, collapse = "|"), names(T21_ge))]
T21_jak <- T21_jak[,1:102]


gtex_jak <- gtex_expression[, grepl(paste(jak_stat$ensembl_gene, collapse = "|"), names(gtex_expression))]
gtex_jak <- gtex_jak[, 1:102]
names(gtex_jak) <- gsub("BLOOD.", "", names(gtex_jak))


write.csv(T21_jak, "recount3_T21.csv")
write.csv(gtex_jak, "recount3_gtex_Blood.csv")


