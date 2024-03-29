#### This program preps the HTP data for downstream analysis
## Pipeline:
##    - load count data
##    - log transform
##    - apply variance filter
##    - pivot to wider

library(tidyverse)
library(caret)
library(sva)
## load expression data
expression <- read.delim("../data/HTP_WholeBlood_RNAseq_FPKMs_Synapse.txt")

## find the number of subjects
n <- length(unique(expression$LabID))

## filter for those with less than 20% missingness
missingness <- expression %>%
    group_by(EnsemblID) %>%
    summarise(missingness = sum(Value ==0)/n) 

miss_filtered <- expression %>%
    filter(EnsemblID %in% (missingness %>% filter(missingness < .2) %>% .$EnsemblID))

## logtranformation
expression$logValue <- log(expression$Value+1)


## filtering based on variance threshold
var_filtered <- miss_filtered %>%
    group_by(EnsemblID) %>%
    summarise(variance = var(logValue)) %>%
    filter(variance > .2)
    

variance_filtered <- miss_filtered %>%
    filter(EnsemblID %in% var_filtered$EnsemblID) %>%
    as_tibble() 


### protein_filtered
prot_filtered <- variance_filtered %>%
    filter(Gene_type == "protein_coding")

## pivot wider. Have to include the unique identifier for row numbers to avoid any errors. 

htp_expr <- prot_filtered %>%
    dplyr::select(LabID, Gene_name,Value) %>%
    distinct(LabID, Gene_name, .keep_all = T) %>%
    pivot_wider(names_from = Gene_name, id_cols = LabID, values_from = Value)




### get the mappings of gene_name, ensembl, and chr from expression
gene_map <- expression %>%
    dplyr::select(Gene_name, EnsemblID,Chr) %>%
    filter(EnsemblID %in% variance_filtered$EnsemblID) %>%
    as_tibble()

expression_list<- list(expression = htp_expr,gene_map = gene_map)
save(expression_list, file = "../data/HTP_transcription_FPKMs_wide_missing_variance_filtered.Rdata")

##### Processing the whole blood transcriptomic data from gtex
library(tidyverse)
library(sva)
load("../data/HTP_transcription_counts_wide_protein_coding_variance_filtered.Rdata")
gtex_expression <- read.csv("../data/gtex/BLOOD_expression.csv")

htp_expr <- expression_list$expression
htp_expr <- as.data.frame(htp_expr)
htp_map <- expression_list$gene_map


htp_map <- htp_map %>%
    distinct(Gene_name, .keep_all = T) %>%
    filter(Gene_name %in% names(htp_expr))

tmp_gtex  <- gtex_expression %>%
    select(c("X", which(gsub("[.].*", "", names(gtex_expression)) %in% gsub("[.].*", "",htp_map$EnsemblID)))) %>%
    column_to_rownames('X')

### GTEX data is RPKMs, divide by 2
tmp_gtex <- as.data.frame(apply(tmp_gtex,2, function(x) x/2))

### Should I add a variance filter to the GTEX data?
var_filtered <- tmp_gtex %>%
    summarise(across(1:ncol(tmp_gtex),  var))

var_filtered <- var_filtered %>%
    select(which(var_filtered > .1))

tmp_gtex <- tmp_gtex %>%
    select(names(var_filtered))

new_names <- htp_map %>%
    filter(gsub("[.].*", "",EnsemblID) %in% gsub("[.].*", "",names(tmp_gtex))) %>%
    arrange(match(gsub("[.].*", "",EnsemblID), gsub("[.].*", "",names(tmp_gtex))))

names(tmp_gtex) <- new_names$Gene_name



### Load in the metadata
htp_meta = read.delim("../../subPhenoDS/data/HTP_data/P4C_metadata_021921_Costello.txt")
htp_meta <- htp_meta %>%
    filter(LabID %in% htp_expr$LabID)
D21_IDs <- htp_meta %>% filter(Karyotype == "Control") %>% .$LabID


### filtering the htp_data
tmp_htp <- htp_expr %>%
    select(c("LabID", which(names(htp_expr) %in% names(tmp_gtex))))%>%
    filter(LabID %in% htp_meta$LabID) %>%
    arrange(match(LabID, htp_meta$LabID)) %>%
    column_to_rownames("LabID") %>%
    relocate(names(tmp_gtex))

### save the sample ids
htp_names <- rownames(tmp_htp)
gtex_names <- rownames(tmp_gtex)

#### Batch correction of the gene expression data
combined <- as.data.frame(t(rbind(tmp_htp, tmp_gtex)))
batch <- as.factor(c(rep(0, nrow(htp_meta)), rep(1, nrow(tmp_gtex))))

library(limma)
combined <- removeBatchEffect(combined, batch)

### combined <- ComBat(dat=combined, batch=batch, mod=NULL, par.prior=TRUE, mean.only=FALSE, ref.batch = 1)

#### scaling the data
##tmp_htp <- apply(tmp_htp,2, scale)
##tmp_gtex <- apply(tmp_gtex, 2,scale )

### prepping for distribution analysis
combined <- as.data.frame(t(combined))
combined <- as.data.frame(apply(combined, 2, scale))
rownames(combined)  <- c(htp_names, gtex_names)

### Separating out the D21 samples 
tmp_combined_D21 <- combined %>%
    filter(rownames(combined) %in% D21_IDs) 

tmp_combined_gtex <- combined %>%
    filter(rownames(combined) %in% gtex_names)

identical(names(tmp_combined_gtex), names(tmp_combined_D21))



### KS test of the data to identify similar distributions

KS <- sapply(1:ncol(tmp_combined_D21), function(x) ks.test(tmp_combined_D21[,x], tmp_combined_gtex[,x])$p.value)
names(KS) <- colnames(tmp_combined_D21)
KS <- sort(KS)
par(mfrow = c(1,1))
plot(KS, main = "K-S pvalues for HTP D21 and GTEX whole blood profiles", xlab = "Gene", pch =16, col = "dark grey")
abline(h = .05, col = "red")





### exploring some of the distributions
# loading the required package
library(dgof)
## s <- names(KS)[length(KS)]
 s <- names(KS)[1]
var1 <- tmp_combined_gtex[,s]
var2 <- tmp_combined_D21[,s]
ks_p <- ks.test(var1, var2)$p.value
# plotting the result
# visualization
par(mfrow = c(2,3))
plot(ecdf(var1),
     xlim = range(c(var1, var2)),
     col = "blue",
     main = paste(s, " - K-S p = ", ks_p ))
plot(ecdf(var2),
     add = TRUE,
     lty = "dashed",
     col = "red")
hist(var1, main = "GTEX")
hist(var2, main = "HTP")
 s <- names(KS)[length(KS)]
var1 <- tmp_combined_gtex[,s]
var2 <- tmp_combined_D21[,s]
ks_p <- ks.test(var1, var2)$p.value
# plotting the result
# visualization
plot(ecdf(var1),
     xlim = range(c(var1, var2)),
     col = "blue",
     main = paste(s, " - K-S p = ", ks_p ))
plot(ecdf(var2),
     add = TRUE,
     lty = "dashed",
     col = "red")
hist(var1, main = "GTEX")
hist(var2, main = "HTP")


### Filter by genes with similar distributions

ks_filt <- names(KS)[KS> .05]
length(ks_filt)

tmp_gtex <- tmp_gtex %>%
    as.tibble() %>%
    dplyr::select(ks_filt)

tmp_htp <- tmp_htp %>%
    as.tibble() %>%
    select(ks_filt)





intrain_htp = createDataPartition(y=htp_meta[,"Karyotype"], p = 0.7)
intrain_gtex = sample(1:nrow(tmp_gtex), .7 * nrow(tmp_gtex))

train_htp <- tmp_htp[intrain$Resample1,]
test_htp <- tmp_htp[-intrain$Resample1,] 

train_gtex <- tmp_gtex[intrain$Resample1,]
test_gtex <- tmp_gtex[-intrain$Resample1,] 

train <- rbind(train_htp, train_gtex)
test <- rbind(test_htp, test_gtex)


y_train <- htp_meta[intrain$Resample1, "Karyotype"]
y_train <- c(y_train, rep("Control", nrow(train_gtex)))

y_test <- htp_meta[-intrain$Resample1, "Karyotype"]
y_test <- c(y_test, rep("Control", nrow(test_gtex)))

y_train <- ifelse(y_train == "T21", 1, 0)
y_test <- ifelse(y_test == "T21", 1, 0)

write.csv(train, "../data/processed/X_train.csv", row.names = F)
write.csv(test, "../data/processed/X_test.csv", row.names = F)
write.csv(y_train, "../data/processed/Y_train.csv", row.names = F)
write.csv(y_test, "../data/processed/Y_test.csv", row.names = F)


