BiocManager::install('limma')


library(limma)
library(openxlsx)

predict11 <- read.csv('./results/T21_blood/predict11.csv', row.names = 1)
predict22 <- read.csv('./results/T21_blood/predict22.csv', row.names = 1)
cross12 <- read.csv('./results/T21_blood/cross12.csv', row.names = 1)
cross21 <- read.csv('./results/T21_blood/cross21.csv', row.names  )

data <- read.csv('./data/recount3_gtex_blood.csv', row.names = 1)
data$type <- 'actual'
names(data) <- names(predict11)

predict11$type = 'predict'
predict22$type = 'predict'
cross12$type = 'cross'
cross21$type = 'cross'


blood_t21 <- rbind(predict22, cross12, data)

mm <- model.matrix(~0 +blood_t21$type)
colnames(mm) <- c("actual","cross", "predict")
fit <- lmFit(t(blood_t21[,1:102]),mm)
contr <- makeContrasts(cross - predict, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit,contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
volcanoplot(tmp, coef = 1, style = "p-value", highlight =20, names = names(fit$Amean), hl.col="blue", xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=.5)

mm <- model.matrix(~0 +blood_t21$type)
colnames(mm) <- c("actual","cross", "predict")
fit <- lmFit(t(blood_t21[,1:102]),mm)
contr <- makeContrasts(actual - predict, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit,contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
volcanoplot(tmp, coef = 1, style = "p-value", highlight =20, names = names(fit$Amean), hl.col="blue", xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=.5)

mm <- model.matrix(~0 +blood_t21$type)
colnames(mm) <- c("actual","cross", "predict")
fit <- lmFit(t(blood_t21[,1:102]),mm)
contr <- makeContrasts(cross-actual, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit,contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

volcanoplot(tmp, coef = 1, style = "p-value", highlight =20, names = names(fit$Amean), hl.col="blue", xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=.5)


fibroblast <- read.xlsx("../../../Grants/INCLUDE/elife-16220-supp1-v2.xlsx", startRow = 2)
fibroblast <- fibroblast[fibroblast$GeneID %in% names(predict11),]
fibroblast <- fibroblast[order(fibroblast$pval),]
fibroblast$rank <- 1:nrow(fibroblast)

exp.top.table <- top.table[rownames(top.table) %in% fibroblast$GeneID,]
exp.top.table$rank <- 1:nrow(exp.top.table)
exp.top.table <- exp.top.table[fibroblast$GeneID,]

cor.test(exp.top.table$rank, fibroblast$rank,method = "kendall")




fibroblast_data <- read.delim("./data/GSE79842_HTSeq_in_fib_counts.txt")


data <- read.csv('./data/recount3_T21.csv', row.names = 1)
data$type <- 'actual'
names(data) <- names(predict11)


blood_t21 <- rbind(predict11, cross21, data)

mm <- model.matrix(~0 +blood_t21$type)
colnames(mm) <- c("actual","cross", "predict")
fit <- lmFit(t(blood_t21[,1:102]),mm)
contr <- makeContrasts(cross - predict, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit,contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
volcanoplot(tmp, coef = 1, style = "p-value", highlight =20, names = names(fit$Amean), hl.col="blue", xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=.5)


mm <- model.matrix(~0 +blood_t21$type)
colnames(mm) <- c("actual","cross", "predict")
fit <- lmFit(t(blood_t21[,1:102]),mm)
contr <- makeContrasts(predict-actual, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit,contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
volcanoplot(tmp, coef = 1, style = "p-value", highlight =20, names = names(fit$Amean), hl.col="blue", xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=.5)

mm <- model.matrix(~0 +blood_t21$type)
colnames(mm) <- c("actual","cross", "predict")
fit <- lmFit(t(blood_t21[,1:102]),mm)
contr <- makeContrasts(cross-actual, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit,contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

volcanoplot(tmp, coef = 1, style = "p-value", highlight =20, names = names(fit$Amean), hl.col="blue", xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=.5)


