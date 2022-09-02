                                        # This is a program for downloading the preliminary data from recount3 into the directory

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("recount3")

## Check that you have a valid Bioconductor installation
BiocManager::valid()

BiocManager::install("LieberInstitute/recount3") # install recount3

library(recount3) # load the recount3 library

human_projects <- available_projects() # load all human projects

dim(human_projects)

head(human_projects)

T21_studies <- c("SRP039348","SRP188973", "SRP188969", "SRP198481", "SRP186520", "SRP186520", "SRP017123", "SRP017123", "SRP017123","SRP017123") # T21 projects


### Potentially follow up here with GTEX samples, the "healthy" ones
