if (!require("gemma.R", quietly = TRUE))
   BiocManager::install("gemma.R")

## load the necessary libraries

library(gemma.R)
library(data.table)
suppressMessages(library(dplyr))

## NB: may be helpful to create a function with a command line argument for efficient downloading of all necessary data

#### Download Down syndrome cohort from GEMMA

## Experiments
DS_experiments <- get_datasets(
   query = "Down Syndrome",
   taxa = "human",
   memoised = getOption("gemma.memoised", TRUE))

save(DS_experiments,file = "./data/GEMMA/DS_experiments.Rdata")

## Gene Expression
DS_expression <- get_dataset_processed_expression(DS_experiments$experiment.ShortName[1],
	      memoised = getOption("gemma.memoised", TRUE))
save(DS_expression, file = "./data/GEMMA/DS_expression.Rdata")

## Sample Metadata
DS_samples <- get_dataset_samples(DS_experiments$experiment.ShortName[1],
	   memoised = getOption("gemma.memoised", TRUE))
save(DS_samples, file = "./data/GEMMA/DS_samples.Rdata")


#### Download AML samples

## Experiments
AML_experiments <- get_datasets(
   query = "GSE15434",
   taxa = "human",
   memoised = getOption("gemma.memoised", TRUE))

save(AML_experiments,file = "./data/GEMMA/AML_experiments.Rdata")

## Gene Expression
AML_expression <- get_dataset_processed_expression(AML_experiments$experiment.ShortName[1],
	      memoised = getOption("gemma.memoised", TRUE))
save(AML_expression, file = "./data/GEMMA/AML_expression.Rdata")

## Sample Metadata
AML_samples <- get_dataset_samples(AML_experiments$experiment.ShortName[1],
	   memoised = getOption("gemma.memoised", TRUE))
save(AML_samples, file = "./data/GEMMA/AML_samples.Rdata")


#### Download PBMC samples

## Experiments
PBMC_experiments <- get_datasets(
   query = "GSE22356",
   taxa = "human",
   memoised = getOption("gemma.memoised", TRUE))

save(PBMC_experiments,file = "./data/GEMMA/PBMC_experiments.Rdata")

## Gene Expression
PBMC_expression <- get_dataset_processed_expression(PBMC_experiments$experiment.ShortName[1],
	      memoised = getOption("gemma.memoised", TRUE))
save(PBMC_expression, file = "./data/GEMMA/PBMC_expression.Rdata")

## Sample Metadata
PBMC_samples <- get_dataset_samples(PBMC_experiments$experiment.ShortName[1],
	   memoised = getOption("gemma.memoised", TRUE))
save(PBMC_samples, file = "./data/GEMMA/PBMC_samples.Rdata")
