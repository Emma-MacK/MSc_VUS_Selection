if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("HDInterval")
install.packages("dplyr")
BiocManager::install("BioQC")
library(BioQC)
library(ggplot2)
library(HDInterval)
library(dplyr)
library(stats)

set.seed(1)  

# want to read in previous functional studies information 
bap1_funt <- readxl::read_excel("extended_data_1_(sge_bap1_dataset).xlsx", sheet = "extended_data_1",  col_names = FALSE,)

# get correct columns name from row 3

colnames(bap1_funt) <- bap1_funt[3,]
bap1_funt <- bap1_funt[-1:-3, ] 


# get rows that coorespond to my variants
 
my_variants <- c("3_52408479_G_A", "3_52409837_C_G", "3_52408046_A_C", 
                 "3_52406252_C_T", "3_52404500_A_C", "3_52409600_C_T")


rows_of_my_variants <- bap1_funt[bap1_funt$chrom_pos_ref_alt %in% my_variants,]

my_variants_mean <-  mean(as.numeric(rows_of_my_variants$functional_score))
# get functional scores only
# get absolute value so it is just change rather than depeleted or enriched

funct_scores <- as.numeric(bap1_funt$functional_score)

# I have 6 variants, so randomly select 6?

# permutation 

n_entries <- length(funct_scores)

bap1_random_selections <- data.frame()

for (x in 1:100000) {
  n <- sample(1:n_entries, 6, replace=FALSE)
  test <- funct_scores[n]
  mean_score <- mean(test)
  mean_diff <- mean(funct_scores) - mean_score
  samples <- append(n,mean_diff)
  bap1_random_selections <- rbind(bap1_random_selections, samples)
}


colnames(bap1_random_selections) <- c("selection 1", "selection 2", 
                                      "selection 3", "selection 4", 
                                      "selection 5", "selection 6",
                                      "mean_diff")


ci95<-hdi(bap1_random_selections$mean_diff)
# what is score of my selection

plot_post_SE1 <- ggplot(bap1_random_selections, aes(x=mean_diff)) +
  geom_histogram(color="darkblue", fill="lightblue", bins = 50) +
  geom_vline(aes(xintercept = mean(bap1_random_selections$mean_diff)), color = "blue") +
  geom_vline(xintercept = as.numeric(ci95[1:2]), linetype="dashed",
             color = "black", size=1, ) +
  geom_vline(aes(xintercept = mean(funct_scores) - my_variants_mean),
             color = "red") +
  ggtitle("Bap 1 random selection mean diff in functional scores")
  
  
plot_post_SE1


# p value 

#  use the distribution to obtain a p-value for our mean-difference by counting how many permuted mean-differences are larger than the one we observed in our actual data. We can then divide this by the number of items in our permutation distribution (i.e., 2000 from our call to replicate, above):
