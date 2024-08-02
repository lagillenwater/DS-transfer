#### This program is for performing pairwise gene evaluation
library(tidyverse)
library("optparse")


option_list = list(
    make_option(c("-f", "--file_1"), type="character", default=NULL, 
              help="First file for pairwise evaluation", metavar="character"),
    make_option(c("-o", "--out_dir"), type="character", default=".", 
                help="output directory [default= current directory]", metavar="character"),
    make_option(c("-s", "--file_2"), type="character", default=NULL, 
                help="second file for pairwise comparison", metavar="character"),
    make_option(c("-s", "--species"), type="character", default='human', 
                help="'human' or 'mouse' [default='human']", metavar="character")    
);



train <- read.csv("../data/processed/X_train.csv")
y_train <- read.csv("../data/processed/Y_train.csv")$x


x_dat <- train[y_train == 1,]
y_dat <- train[y_train == 0,]

## column medians
gene_medians <- x_dat %>%
    summarise_if(is.numeric, mean) %>%
    t() %>%
    as.numeric()

## find gene ratios

gene_ratios <- outer(gene_medians ,gene_medians, "/")
colnames(gene_ratios) <- names(x_dat)
rownames(gene_ratios) <- names(x_dat)

## pivot to longer
test <- gene_ratios %>%
    as_tibble() %>%
    mutate(Gene1 = rownames(gene_ratios))%>%
    pivot_longer(cols = colnames(gene_ratios),  names_to = "Gene2", values_to = "gene_ratio") %>%
    filter(Gene1 != Gene2)



## column medians
gene_medians <- y_dat %>%
    summarise_if(is.numeric, mean) %>%
    t() %>%
    as.numeric()

## find gene ratios

gene_ratios <- outer(gene_medians ,gene_medians, "/")
colnames(gene_ratios) <- names(x_dat)
rownames(gene_ratios) <- names(x_dat)

## pivot to longer
test2 <- gene_ratios %>%
    as_tibble() %>%
    mutate(Gene1 = rownames(gene_ratios))%>%
    pivot_longer(cols = colnames(gene_ratios),  names_to = "Gene2", values_to = "gene_ratio") %>%
    filter(Gene1 != Gene2)


plot(test[1:100,"gene_ratio"]$gene_ratio,test2[1:100,"gene_ratio"]$gene_ratio, xlim = c(-50,50), ylim = c(-50,50))

cor(test[,"gene_ratio"]$gene_ratio,test2[,"gene_ratio"]$gene_ratio)
