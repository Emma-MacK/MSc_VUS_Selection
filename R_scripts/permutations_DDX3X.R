library(BioQC)
library(ggplot2)
library(HDInterval)
library(dplyr)
library(stringr)
library(purrr)


set.seed(1) 

# want to read in previous functional studies information 

file_list <- list.files(path = "DDX3X_functional_scores/", pattern = "*.csv", full.names = TRUE)
data_list <- lapply(file_list, read.csv)

# Optionally, combine all data frames into one
DDX3X_funct <- do.call(rbind, data_list)

my_variants <- c("c.1451T>C", "c.1487_1488del", "c.931C>T",
                 "c.1126C>T", "c.611del", "c.1541T>C", "c.1724C>T",
                 "c.1660A>T", "c.1301T>A", "c.1126C>T", "c.603_604del",
                 "c.1314A>G", "c.691dup", "c.1582C>T", "c.1598_1601dup",
                 "c.769A>C", "c.824C>T", "c.1236dup", "c.533A>G",
                 "c.1457_1458insT", "c.441del", "c.455_456dup", "c.1315+1G>A",
                 "c.1490C>T", "c.1490C>A", "c.1382T>C", "c.1600C>T",
                 "c.1865G>A", "c.1445_1447del", "c.295C>T", "c.1052G>A",
                 "c.146C>T", "c.1319A>G", "c.536T>C", "c.146C>T")

# get rows that coorespond to my variants



hgvs <- DDX3X_funct$hgvs_splice
my_variants_IDs <- hgvs[str_detect(hgvs, paste(my_variants, collapse = "|"))]
my_variants_IDs <- na.omit(my_variants_IDs)

my_variant_data <- DDX3X_funct[DDX3X_funct$hgvs_splice %in% my_variants_IDs,]

my_variants_mean <-  mean(as.numeric(my_variant_data$d2))
# only 27 as 146 appears twice

# get functional scores only


#  Variants with cLFC FDR>0.01 on Day 15 were classified as "SGE-unchanged".



funct_scores <- as.numeric(DDX3X_funct$D21_combined_LFC)


funct_scores_plot <- ggplot(DDX3X_funct, aes(x=as.numeric(DDX3X_funct$D21_combined_LFC))) +
  geom_histogram(color="darkblue", fill="lightblue", bins = 50) +
  geom_boxplot(data = DDX3X_funct[DDX3X_funct$hgvs_splice %in% my_variants_IDs,],
               aes(x = as.numeric(my_variant_data$D21_combined_LFC)), 
               fill = "grey", alpha = 75, width = 100, position = position_nudge(y = 5000)) +
  geom_vline(aes(xintercept = mean(as.numeric(DDX3X_funct$D21_combined_LFC))), color = "black") +
  ggtitle("combined LFC of DDX3X variants at day 21") +
  xlab("combined LFC at day 21") +
  theme_bw()

funct_scores_plot

fuct_scores_FDR <- ggplot(DDX3X_funct, aes(x=as.numeric(D21_combined_LFC), y = log(as.numeric(cLFCd21_BH_FDR)))) +
  geom_point(color="lightblue", fill="lightblue") +
  geom_point(data = DDX3X_funct[DDX3X_funct$hgvs_splice %in% my_variants_IDs,], color="red")+
  ggtitle("DDX3X Day 21 functional score v FDR") +
  xlab("functional score (Log FC) at Day 21") +
  ylab("log 10 FDR") +
  scale_y_reverse() +
  theme_bw() +
  geom_hline(aes(yintercept = log(0.01)), color = "darkgrey") 



fuct_scores_FDR
# I have 6 variants, so randomly select 6?

# permutation 

n_entries <- length(funct_scores)

DDX3X_random_selections <- data.frame()

for (x in 1:10000) {
  n <- sample(1:n_entries, 27, replace=FALSE)
  test <- funct_scores[n]
  mean_score <- mean(test)
  samples <- append(n,mean_score)
  DDX3X_random_selections <- rbind(DDX3X_random_selections, samples)
}

colnames(DDX3X_random_selections) <- c("selection 1", "selection 2", 
                                      "selection 3", "selection 4", 
                                      "selection 5", "selection 6",
                                      "selection 7", "selection 8",
                                      "selection 9", "selection 10",
                                      "selection 11", "selection 12", 
                                      "selection 13", "selection 14", 
                                      "selection 15", "selection 16",
                                      "selection 17", "selection 18",
                                      "selection 19", "selection 20",
                                      "selection 21", "selection 22", 
                                      "selection 23", "selection 24", 
                                      "selection 25", "selection 26",
                                      "selection 27",
                                      "mean_score")



pos_3sd = mean(DDX3X_random_selections$mean_score) + 3*sd(DDX3X_random_selections$mean_score)
neg_3sd = mean(DDX3X_random_selections$mean_score) - 3*sd(DDX3X_random_selections$mean_score)
# what is score of my selection

plot_post_SE1 <- ggplot(DDX3X_random_selections, aes(x=mean_score)) +
  geom_histogram(color="darkblue", fill="lightblue", bins = 30) +
  geom_vline(aes(xintercept = mean(DDX3X_random_selections$mean_score)), color = "black") +
  geom_vline(xintercept =c(pos_3sd, neg_3sd), linetype="dashed",
             color = "black", size=1, ) +
  xlab("mean cLFC") +
  geom_vline(aes(xintercept = mean(my_variant_data$D21_combined_LFC)), color = "red") +
  ggtitle("DDX3X - mean cLFC at day 21 of random selection ") +
  theme_bw()


plot_post_SE1

std.error <- function(x) sd(x)/sqrt(length(x))

pop_mean <- mean(as.numeric(DDX3X_funct$D21_combined_LFC))
t.test(as.numeric(my_variant_data$D21_combined_LFC), mu= pop_mean)

t.test(as.numeric(my_variant_data$D21_combined_LFC),DDX3X_random_selections$mean_score)
