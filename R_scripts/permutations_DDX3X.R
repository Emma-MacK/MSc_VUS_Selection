library(BioQC)
library(ggplot2)
library(HDInterval)
library(dplyr)
library(stringr)
library(purrr)

set.seed(1) 

# want to read in previous functional studies information 
# DDX3X dataset has functional scores for each exon in seperate file

file_list <- list.files(path = "DDX3X_functional_scores/", pattern = "*.csv", full.names = TRUE)
data_list <- lapply(file_list, read.csv)

# combine all data frames into one
DDX3X_funt <- do.call(rbind, data_list)

# get rows that coorespond to my variants - have previously checked transcripts are comparible and same protein change
my_variants <- c("c.1724C>T", "c.1301T>A", "c.1314A>G", 
                 "c.533A>G", "c.1865G>A", "c.295C>T",
                 "c.146C>T", "c.1319A>G", "c.536T>C",
                 "c.146C>T")

hgvs <- DDX3X_funt$hgvs_splice
my_variants_IDs <- hgvs[str_detect(hgvs, paste(my_variants, collapse = "|"))]
# selecting some NA columns
my_variants_IDs <- na.omit(my_variants_IDs)
my_variant_data <- DDX3X_funt[DDX3X_funt$hgvs_splice %in% my_variants_IDs,]

my_variants_mean <-  mean(as.numeric(my_variant_data$score))
# only 9 as 146 appears twice 

# get functional scores only
# get absolute value so it is just change rather than depeleted or enriched

funct_scores <- as.numeric(DDX3X_funt$score)

# permutation 
n_entries <- length(funct_scores)
DDX3X_random_selections <- data.frame()
for (x in 1:100000) {
  n <- sample(1:n_entries, 10, replace=FALSE)
  test <- funct_scores[n]
  mean_score <- mean(test)
  mean_diff <- mean(funct_scores) - mean_score
  samples <- append(n,mean_diff)
  DDX3X_random_selections <- rbind(DDX3X_random_selections, samples)
}

colnames(DDX3X_random_selections) <- c("selection 1", "selection 2", 
                                      "selection 3", "selection 4", 
                                      "selection 5", "selection 6",
                                      "selection 7", "selection 8",
                                      "selection 9", "selection 10",
                                      "mean_diff")

# calculate confidence intervals 
ci95<-hdi(DDX3X_random_selections$mean)

# plot 
plot_DDX3X_mean_diff <- ggplot(DDX3X_random_selections, aes(x=mean_diff)) +
  geom_histogram(color="darkblue", fill="lightblue", bins = 50) +
  geom_vline(aes(xintercept = mean(DDX3X_random_selections$mean_diff)), color = "blue") +
  geom_vline(xintercept = as.numeric(ci95[1:2]), linetype="dashed",
             color = "black", size=1, ) +
  geom_vline(aes(xintercept = mean(funct_scores) - my_variants_mean), color = "red") +
  ggtitle("DDX3X mean diff permutations")

plot_DDX3X_mean_diff



