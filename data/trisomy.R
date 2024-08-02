# R script to download selected samples
# Copy code and run on a local machine to initiate download

library("rhdf5")    # can be installed using Bioconductor
options(timeout = 1e9)
destination_file = "human_gene_v2.2.h5"
extracted_expression_file = "trisomy_expression_matrix.tsv"
url = "https://s3.dev.maayanlab.cloud/archs4/files/human_gene_v2.2.h5"

# Check if gene expression file was already downloaded, if not in current directory download file form repository
if(!file.exists(destination_file)){
    print("Downloading compressed gene expression matrix.")
    download.file(url, destination_file, quiet = FALSE, mode = 'wb')
}

# Selected samples to be extracted
samp = c("GSM1261907","GSM1261908","GSM1261906","GSM2676183","GSM3144757","GSM3144758","GSM3144759","GSM3144760","GSM3144761","GSM3144762","GSM3144763","GSM3144764","GSM3144765","GSM3144766","GSM3144767","GSM3144768","GSM3144769","GSM3144770","GSM3144771","GSM3144772","GSM3144773","GSM3144774","GSM3144775","GSM3144776","GSM3144777","GSM3144778","GSM3144779","GSM3144780","GSM3144781","GSM3144782","GSM3144783",
"GSM3144784","GSM3144785","GSM3144786","GSM3144787","GSM3144788","GSM3144789","GSM3144790","GSM3144791","GSM3144792","GSM3144793","GSM3144794","GSM3144795","GSM3144796","GSM3144797","GSM3144798","GSM3144799","GSM3144800","GSM3144801","GSM3144802","GSM3144803","GSM3144804","GSM3144805","GSM3144806","GSM3144807","GSM3144808","GSM3144809","")

# Retrieve information from compressed data
samples = h5read(destination_file, "meta/samples/geo_accession")
genes = h5read(destination_file, "meta/genes/symbol")

# Identify columns to be extracted
sample_locations = which(samples %in% samp)

# extract gene expression from compressed data
expression = t(h5read(destination_file, "data/expression", index=list(sample_locations, 1:length(genes))))
rownames(expression) = genes
colnames(expression) = samples[sample_locations]


rownames(expression) = genes
colnames(expression) = samples[sample_locations]

H5close()


# Print file
write.table(expression, file=extracted_expression_file, sep="\t", quote=FALSE, col.names=NA)
print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_file))

## extract geo metadata
library(GEOquery)

meta <- lapply(c(sample_locations, sample(1:length(samples), 100)), function(x) Meta(getGEO(samples[x]))[c("geo_accession","characteristics_ch1", "source_name_ch1")])

# write unstructured metadata to file
sink("metadata.txt")
print(meta)
sink()

