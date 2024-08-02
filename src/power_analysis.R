install.packages('pwr') # install pwr package
BiocManager::install('RnaSeqSampleSize')

library(pwr)
library(tidyverse)
library(ggplot2)
library(RnaSeqSampleSize)

comorbidities <- read.delim("./data/P4C_Comorbidity_020921.tsv", skip = 1)
head(comorbidities)

metadata <-read.delim("./data/P4C_metadata_021921_Costello.txt")
head(metadata)

comorbidities <- comorbidities %>%
    select(RecordID, Condition, HasCondition) %>%
    pivot_wider(names_from = Condition, values_from = HasCondition)


meta <- metadata %>%
    inner_join(comorbidities, by = "RecordID")
head(meta)

ds_meta <- meta %>%
    filter(Karyotype == "T21")
head(ds_meta)
dim(ds_meta) # 384 subjects with DS

sapply(names(comorbidities)[2:ncol(comorbidities)], function(x) table( ds_meta[[x]]))


                                        # power analysis for case imbalanced data


est_power_curve(n = 100, w = 1, lambda0 = 1)

powers <- seq(.1,.9,.2)
ratios <- seq(.1,.5,.1)


sample_power <- matrix(0,ncol = 6, nrow = 5)
sample_power[,1] = ratios

for(power in 1:length(powers)){
    print(power)
    sample_size_vector <- numeric()
    for( k in ratios ){
        print(k)
        sample_size_vector<-  c(sample_size_vector,sample_size(power = powers[power],m1 = 164, m = 16411,rho=2,  k = k, f = .05, phi0 = 1,lambda0 =5 ))
    }
    sample_power[,1+power] <- sample_size_vector
}
sample_power <- as.data.frame(sample_power)
names(sample_power) <- c("ratio",  powers)
sample_power <- sample_power %>%
    pivot_longer(!ratio, names_to = "power", values_to = "sample_size") 
sample_power$power <- as.numeric(sample_power$power)

# plot of power curves
p1 <- ggplot(sample_power, aes(x = sample_size, y = power, color = factor(ratio))) +
    geom_line(size = 2) +
    geom_point(size = 3) +
    ylim(c(0,1)) +
    theme_classic()+
    geom_hline(yintercept = .8, linetype = "dotted")+
    xlab("sample size") +
    guides(color = guide_legend(title = "Ratio of \nCases/Controls")) +
    theme(legend.position = c(0.9,.5))

pdf("./results/figures/power_analysis.pdf")
print(p1)
dev.off()
                     

load("./data/DS_transcription_profiles_wide.RData")

ds_wide <- ds_wide %>%
    distinct(LabID, .keep_all = T) # remove duplicate names


                                        # separate by karyotype

ds_wide_ds <- ds_wide %>%
    filter(LabID %in% ds_meta$LabID)
#variance filtere
ds_variance <- apply(ds_wide_ds[,2:ncol(ds_wide_ds)],2,var)

                                        # numebr of genes in DS cohort with VAR > .1

ds_wide_var_filt <- ds_wide_ds %>%
    select(names(ds_variance[(ds_variance > .1)]))



                                        # Identify subsets of disease
ds_meta <- ds_meta %>%
    filter(LabID %in% ds_wide_ds$LabID)

names(ds_meta) <- make.names(names(ds_meta))

ds_meta %>%
    filter(Any.autoimmune.skin.condition == 1 ) # 123

ds_meta %>%
    filter(Any.autoimmune.skin.condition == 1 & Any.sleep.apnea == 1  ) # 74

ds_meta %>%
    filter(Any.autoimmune.skin.condition == 1 & Any.sleep.apnea == 1 & Obesity == 1 ) # 38


sample_size(power = .8,m1 = 164, m = 16411,rho=2,  k = .42, f = .05, phi0 = 1,lambda0 =5 )
sample_size(power = .8,m1 = 164, m = 16411,rho=2,  k = .26, f = .05, phi0 = 1,lambda0 =5 )
sample_size(power = .8,m1 = 164, m = 16411,rho=2,  k = .13, f = .05, phi0 = 1,lambda0 =5 )



### exploration of intestine-specific network
intestine_graph <- read.delim("./data/intestine_top", header = F)




####Drug Bank
biokg <- read.delim("./biokg/biokg.links.tsv", header =F)
protein <- read.delim("./biokg/biokg.metadata.protein.tsv", header = F)
drug_meta <- read.delim("./biokg/biokg.metadata.drug.tsv", header = F)
head(drug_meta)
table(drug_meta$V2)
mtx <- drug_meta %>%
    filter(V3 == "Methotrexate")
mtx <- "DB00563"


mtx_graph <- biokg %>%
    filter(V1 == mtx) %>%
    filter(!(V2 %in%  c("DPI", "DDI", "DRUG_DISEASE_ASSOCIATION", "DRUG_PATHWAY_ASSOCIATION")))
noquote(mtx_graph$V3)



    


dpz <- drug_meta %>%
    filter(V3 == "Donepezil")
dpz <- "DB00843"

dpz_graph <- biokg %>%
    filter(V1 == dpz) %>%
    filter(!(V2 %in%  c("DPI", "DDI", "DRUG_DISEASE_ASSOCIATION", "DRUG_PATHWAY_ASSOCIATION"))) 
dpz_graph

dpz_genes <- dpz_graph %>%
    .$V3 %>%
    noquote


library(openxlsx)
deg <- read.xlsx("./data/elife-16220-supp1-v2.xlsx", sheet =2, startRow=2) # lymphoblast cell lines


dpz_genes <- read.delim("./data/dpz_genes")

dpz_genes <- dpz_genes %>%
    inner_join(deg, by = c("To" = "GeneID"))

dpz_genes %>% inner_join(dpz_graph, by = c("From" = "V3"))

                                        # List of proteins involved



                                        # intestitnal regulatory network
library(igraph)
intestine <- read.delim("./data/intestine_top", header = F)
int_graph <- graph_from_data_frame(intestine,directed = F)


    
BiocManager::install('UniProt.ws')
BiocManager::install("org.Hs.eg.db")
BiocManager::install('RCy3')
library(UniProt.ws)
library(org.Hs.eg.db)
library(RCy3)

mtx_entrez_ids <- select(org.Hs.eg.db, mtx_graph$V3, "ENTREZID", "UNIPROT")
mtx_graph <-  mtx_graph %>%
    inner_join(mtx_entrez_ids, by = c("V3" = "UNIPROT"))


mtx_transporters <- mtx_graph %>%
    filter(V2 == "DRUG_TRANSPORTER")

mtx_transporter_neighbors <- mtx_transporters$ENTREZID
for(gene in mtx_transporters$ENTREZID){
    mtx_transporter_neighbors <- c(mtx_transporter_neighbors, names(neighbors(int_graph, gene)))
    }
    
mtx_transporter_int_graph <- subgraph(int_graph, mtx_transporter_neighbors)
summary(mtx_transporter_int_graph)

createNetworkFromIgraph(mtx_transporter_int_graph)

load("./data/GTEx_PANDA_tissues.RData")

head(genes)
genes %>% filter(Symbol == "SLC19A1")


