##### This program is for an analysis comparing the distributions between two separate GTEX datasets downloaded from the recount3 database.
##### First comparing blood and skin
library(tidyverse)
library(sva)
blood <- read.csv("../data/gtex/BLOOD_expression.csv")

skin <- read.csv("../data/gtex/skin_expression.csv")

### check for identical names
identical(names(blood), names(skin))

### variance filter
blood_var <- blood  %>%
    summarise(across(2:ncol(blood),  var))

skin_var <- skin  %>%
    summarise(across(2:ncol(skin),  var))


tokeep <- union(names(blood_var)[blood_var > .1],names(skin_var)[skin_var > .1])


tmp_blood <- blood %>%
    select(tokeep)

tmp_skin <- skin %>%
    select(tokeep)



#### Compare distributions prior to scaling
KS <- sapply(1:ncol(tmp_blood), function(x) ks.test(tmp_blood[,x], tmp_skin[,x])$p.value)
names(KS) <- colnames(tmp_blood)
KS <- sort(KS)
par(mfrow = c(1,1))
plot(KS, main = "K-S pvalues for GTEX blood and skin profiles", xlab = "Gene", pch =16, col = "dark grey")
abline(h = .05, col = "red")



### exploring some of the distributions
# loading the required package
library(dgof)
## s <- names(KS)[length(KS)]
 s <- names(KS)[1]
var1 <-tmp_blood[,s]
var2 <- tmp_skin[,s]
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
hist(var1, main = "blood")
hist(var2, main = "skin")
 s <- names(KS)[length(KS)-100]
var1 <- tmp_blood[,s]
var2 <- tmp_skin[,s]
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
hist(var1, main = "blood")
hist(var2, main = "skin")

### save the sample ids
htp_names <- rownames(tmp_htp)
gtex_names <- rownames(tmp_gtex)

#### Batch correction of the gene expression data
combined <- as.data.frame(t(rbind(tmp_blood, tmp_skin)))
batch <- c(rep(0, nrow(tmp_blood)), rep(1, nrow(tmp_skin)))
combined <- ComBat(dat=combined, batch=batch, mod=NULL, par.prior=TRUE, mean.only=FALSE, ref.batch = 0)

#### scaling the data
##tmp_htp <- apply(tmp_htp,2, scale)
##tmp_gtex <- apply(tmp_gtex, 2,scale )

### prepping for distribution analysis
combined <- as.data.frame(t(combined))
combined <- as.data.frame(apply(combined, 2, scale))
rownames(combined)  <- c(blood$X, skin$X)

### Separating out the D21 samples 
tmp_combined_blood <- combined %>%
    filter(rownames(combined) %in% blood$X) 



tmp_combined_skin <- combined %>%
    filter(rownames(combined) %in% skin$X)


#### Compare distributions prior to scaling
KS <- sapply(1:ncol(tmp_combined_blood), function(x) ks.test(tmp_combined_blood[,x], tmp_combined_skin[,x])$p.value)
names(KS) <- colnames(tmp_combined_blood)
KS <- sort(KS)
par(mfrow = c(1,1))
plot(KS, main = "K-S pvalues for GTEX blood and skin profiles", xlab = "Gene", pch =16, col = "dark grey")
abline(h = .05, col = "red")



### exploring some of the distributions
# loading the required package
library(dgof)
## s <- names(KS)[length(KS)]
 s <- names(KS)[1]
var1 <-tmp_combined_blood[,s]
var2 <- tmp_combined_skin[,s]
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
hist(var1, main = "blood")
hist(var2, main = "skin")
 s <- names(KS)[length(KS)-3]
var1 <- tmp_combined_blood[,s]
var2 <- tmp_combined_skin[,s]
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
hist(var1, main = "blood")
hist(var2, main = "skin")


