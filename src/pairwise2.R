library(tidyverse)
library(ggplot2)


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


load("./data/DS_transcription_profiles_wide.RData")

ds_wide <- ds_wide %>%
    distinct(LabID, .keep_all = T) # remove duplicate names


                                        # separate by karyotype

ds_wide_T21 <- ds_wide %>%
    filter(LabID %in% ds_meta$LabID)

ds_wide_D21 <- ds_wide %>%
    filter(!(LabID %in% ds_meta$LabID))


## variance filter function
varFilter <- function(gene_expression,threshold) {
    #variance filtere
    variance <- apply(gene_expression[,2:ncol(gene_expression)],2,var)
                                        # numebr of genes in DS cohort with VAR > .1
    filtered <- gene_expression %>%
        dplyr::select(LabID,names(variance[(variance > threshold)]))
    return(filtered)
}

ds_T21 <- varFilter(ds_wide_T21, .1)
ds_D21 <- varFilter(ds_wide_D21, .1)

ds_T21 <- ds_T21[, names(ds_T21) %in% names(ds_D21)]
ds_D21 <- ds_D21[, names(ds_D21) %in% names(ds_T21)]

ds_T21[,2:ncol(ds_T21)] <- apply(ds_T21[,2:ncol(ds_T21)], 2, log)
ds_D21[,2:ncol(ds_D21)] <- apply(ds_D21[,2:ncol(ds_D21)], 2, log)

library(openxlsx)
chr21 <- read.xlsx("./data/CHR21.xlsx")
chr21 <- chr21 %>% filter(e.symbol %in% names(ds_T21))
genes <- c(which(names(ds_T21) %in% chr21$e.symbol), sample(2:ncol(ds_T21), 59))
genes <- unique(genes)

pairwiseRatios <- function(gene_expression, genes) {

    gene_expression <- gene_expression[,genes]
    ratios <- sapply(2:(ncol(gene_expression)-1), function(i) {
        print(i/ncol(gene_expression-1))
        sapply((i+1):ncol(gene_expression), function(j) {
            x = median(unlist(gene_expression[,i]/gene_expression[,j]))
            names(x) = paste(names(gene_expression)[i], names(gene_expression)[j])
            return(x)
        }, USE.NAMES = TRUE)})
    ratios <- unlist(ratios)
     return(ratios)
}

T21_ratios <- pairwiseRatios(ds_T21,genes)
D21_ratios <- pairwiseRatios(ds_D21, genes)





T21_ratios_noinf <- T21_ratios[names(T21_ratios) %in% intersect(names(T21_ratios)[is.finite(T21_ratios)] , names(D21_ratios)[is.finite(D21_ratios)])]
D21_ratios_noinf <- D21_ratios[names(D21_ratios) %in% intersect(names(T21_ratios)[is.finite(T21_ratios)] , names(D21_ratios)[is.finite(D21_ratios)])]



head(T21_ratios_noinf[order(T21_ratios_noinf, decreasing = T)])
head(D21_ratios_noinf[order(D21_ratios_noinf, decreasing = T)])

toplot <- as.data.frame(cbind(T21_ratios_noinf, D21_ratios_noinf))

fun <- function(x, y) {
    grepl(x, y)
}



toplot$yes_chr21 <- ifelse(mapply(fun,paste(str_split(rownames(toplot), " ")), chr21$e.symbol), "yes","no")



ggplot(data = toplot,  aes(x = T21_ratios_noinf, y = D21_ratios_noinf)) +
    geom_point(size = 3)+
    theme_classic() +
    geom_abline(intercept = 0, linetype='dashed', col = 'red', size = 2) +
    xlim(c(-10,10))+
    ylim(c(-10,10)) 
    

rand <- rnorm(nrow(toplot), mean = 0, sd = .3)
toplot$rand <- toplot$T21_ratios_noinf+ rand


ggplot(data = toplot,  aes(x = T21_ratios_noinf, y = rand)) +
    geom_point(size = 3)+
    theme_classic()+
    geom_abline(intercept = 0, linetype='dashed', col = 'red', size = 2) +
    xlim(c(-10,10))+
    ylim(c(-10,10)) 
    


names(meta) <- make.names(names(meta))

T21_skin <- ds_T21 %>%
    filter(LabID %in% (meta %>% filter(Any.autoimmune.skin.condition ==1) %>% .$LabID))
D21_noskin <- ds_T21 %>%
    filter(LabID %in% (meta %>% filter(!(Any.autoimmune.skin.condition ==1)) %>% .$LabID))


T21_skin_ratios <- pairwiseRatios(T21_skin, genes)
D21_skin_ratios <- pairwiseRatios(D21_skin, genes)
T21_skin_ratios_noinf <- T21_skin_ratios[is.finite(T21_skin_ratios) & is.finite(D21_skin_ratios)]
D21_skin_ratios_noinf <- D21_skin_ratios[is.finite(T21_skin_ratios) & is.finite(D21_skin_ratios)]


T21_ratios_noinf <- T21_ratios[names(T21_ratios) %in% union(names(T21_ratios)[is.finite(T21_ratios)] , names(T21_skin_ratios)[is.finite(T21_skin_ratios)])]
T21_skin_ratios_noinf <- T21_skin_ratios[names(T21_skin_ratios) %in% union(names(T21_ratios)[is.finite(T21_ratios)] , names(T21_skin_ratios)[is.finite(T21_skin_ratios)])]
T21_skin_ratios_noinf <- T21_skin_ratios_noinf[names(T21_ratios_noinf)]


toplot2 <- as.data.frame(cbind(T21_ratios_noinf, T21_skin_ratios_noinf))
ggplot(data = toplot2,  aes(x = T21_ratios_noinf, y = T21_skin_ratios_noinf)) +
    geom_point() +
    theme_classic()

rand <- rnorm(nrow(toplot2), mean = .05, sd = .00001)
toplot2$rand <- T21_ratios+rand
ggplot(data = toplot2,  aes(x = T21_ratios, y = rand)) +
    geom_point() +
    theme_classic()




library(openxlsx)
library('RColorBrewer')

deg <- read.xlsx("./data/elife-16220-supp1-v2.xlsx", sheet =2, startRow=2) # lymphoblast cell lines
toplot3 <- head(deg, 20)
ggplot(toplot3, aes(x = factor(1),y = GeneID,fill = log2FoldChange_adj )) +
    geom_tile() +
    scale_fill_gradient2(midpoint = mean(toplot3$log2FoldChange_adj), low = "blue", mid = "white", high = "red")+
    theme_minimal()
    



