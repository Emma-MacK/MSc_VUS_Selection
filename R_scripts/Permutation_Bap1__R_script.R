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
                 "3_52408062_A_G", "3_52408061_C_G",
                 "3_52406252_C_T", "3_52404500_A_C", "3_52409600_C_T",
                 "3_52408479_G_A")


rows_of_my_variants <- bap1_funt[bap1_funt$chrom_pos_ref_alt %in% my_variants,]

my_variants_mean <-  mean(as.numeric(rows_of_my_variants$functional_score))
# get functional scores only
# get absolute value so it is just change rather than depeleted or enriched

funct_scores <- as.numeric(bap1_funt$processed_LFC_D4_D21)
pop_mean <- mean(funct_scores)
# I have 6 variants, so randomly select 6?

pop_mean <- mean(as.numeric(bap1_funt$processed_LFC_D4_D21))
t.test(as.numeric(bap1_funt$processed_LFC_D4_D21), as.numeric(rows_of_my_variants$processed_LFC_D4_D21) )

#bap1_funt$processed_BH_FDR_D4_D21
#bap1_funt$processed_LFC_D4_D21


funct_scores_plot <- ggplot(bap1_funt, aes(x=as.numeric(processed_LFC_D4_D21))) +
  geom_histogram(color="darkblue", fill="lightblue", bins = 50) +
  geom_boxplot(data = bap1_funt[bap1_funt$chrom_pos_ref_alt %in% my_variants,],
               aes(x = as.numeric(rows_of_my_variants$processed_LFC_D4_D21)), 
               fill = "grey", alpha = 75, width = 100, position = position_nudge(y = 5000)) +
  geom_vline(aes(xintercept = mean(as.numeric(processed_LFC_D4_D21))), color = "black") +
  ggtitle("combined LFC of BAP1 variants at day 21") +
  xlab("combined LFC at day 21") +
  theme_bw()

funct_scores_plot

fuct_scores_FDR <- ggplot(bap1_funt, aes(x=as.numeric(processed_LFC_D4_D21), y = log(as.numeric(processed_BH_FDR_D4_D21)))) +
  geom_point(color="lightblue", fill="lightblue") +
  geom_point(data = bap1_funt[bap1_funt$chrom_pos_ref_alt %in% my_variants,], color="red")+
  ggtitle("BAP1 - Day 21 functional score v FDR") +
  xlab("functional score (Log FC) at Day 21") +
  ylab("log FDR") +
  scale_y_reverse() +
  theme_bw() +
  geom_hline(aes(yintercept = log(0.01)), color = "darkgrey")   

fuct_scores_FDR


# permutation 

n_entries <- length(funct_scores)

bap1_random_selections <- data.frame()

for (x in 1:10000) {
  n <- sample(1:n_entries, 8, replace=FALSE)
  test <- funct_scores[n]
  mean_score <- mean(test)
  samples <- append(n,mean_score)
  bap1_random_selections <- rbind(bap1_random_selections, samples)
}



colnames(bap1_random_selections) <- c("selection 1", "selection 2", 
                                      "selection 3", "selection 4", 
                                      "selection 5", "selection 6",
                                      "selection 7", "selection 8", "mean_score")


pos_3sd = mean(bap1_random_selections$mean_score) + 3*sd(bap1_random_selections$mean_score)
neg_3sd = mean(bap1_random_selections$mean_score) - 3*sd(bap1_random_selections$mean_score)
# what is score of my selection

plot_post_SE1 <- ggplot(bap1_random_selections, aes(x=mean_score)) +
  geom_histogram(color="darkblue", fill="lightblue", bins = 30) +
  geom_vline(aes(xintercept = mean(as.numeric(bap1_random_selections$mean_score))), color = "black") +
  geom_vline(xintercept = c(pos_3sd,neg_3sd), linetype="dashed",
             color = "black", size=1, ) +
  xlab("mean cLFC") +
  geom_vline(aes(xintercept = mean(as.numeric(rows_of_my_variants$processed_LFC_D4_D21))),
             color = "red") +
  ggtitle("BAP1 - mean cLFC at day 21 of random selections") +
  theme_bw()
  
  
plot_post_SE1
# t.test(as.numeric(rows_of_my_variants$functional_score), bap1_random_selections$mean_score)

# p value 

# try without known pathogenic
# !completecases needed as rows with no clinvar entry returned as all NA otherwise
BAP_funt_no_path <- bap1_funt[!bap1_funt$clinvar_clinical_significance == "Pathogenic"  & !bap1_funt$clinvar_clinical_significance == "Likely pathogenic"| !complete.cases(bap1_funt$clinvar_clinical_significance),]
my_variants_no_path <- rows_of_my_variants[!rows_of_my_variants$clinvar_clinical_significance =="Pathogenic" & !rows_of_my_variants$clinvar_clinical_significance == "Likely pathogenic"| !complete.cases((rows_of_my_variants$clinvar_clinical_significance)),] 

mean_no_path <- mean(as.numeric(BAP_funt_no_path$functional_score))
my_var_no_path <- as.numeric(my_variants_no_path$functional_score)




#  use the distribution to obtain a p-value for our mean-difference by counting how many permuted mean-differences are larger than the one we observed in our actual data. We can then divide this by the number of items in our permutation distribution (i.e., 2000 from our call to replicate, above):

