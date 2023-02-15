### 2023-02-09
### Program for downloading the data file from cavatica.
### Created a list of download links from Cavatica. Files had been linked from INCLUDE
setwd("./data/HTP_raw_counts")
files <- read.delim("download-links.txt", header = F)
command <- lapply(files, function(x) paste0("wget -O ", gsub(".*[/]", "",gsub('[?].*', '', x))," '", x, "'"))
sapply(command[[1]], function(x) system(x))
