
import sys
import numpy as np
import functools
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import vcf

# get gene name
# read in vcf and csv outputs

exomiser_csv=pd.read_csv("all_genes_merged_variants.csv")
filtered_csv=pd.read_csv("all_genes_filtered_variants_stage_4.csv")
# violin plot of exomiser score different variant types

# remove columns without numeric score
exomiser_csv['score'] = pd.to_numeric(exomiser_csv['score'], errors='coerce')
filtered_csv['score'] = pd.to_numeric(filtered_csv['score'], errors='coerce')

participants = exomiser_csv.iloc[:, [1]]

participants.to_csv("participants_stage_4.csv")

#print(exomiser_csv['score'])
#print(exomiser_csv)
#
exomiser_csv = exomiser_csv.dropna(subset=['score'])
filtered_csv = filtered_csv.dropna(subset=['score'])


#print(exomiser_csv['score'])
#print(exomiser_csv)

# add NAN where clivar is not available
no_class_str = "No classification"

exomiser_csv['clinvar'].fillna(no_class_str, inplace=True)
filtered_csv['clinvar'].fillna(no_class_str, inplace=True)


# set up dataframe 

fig, axs = plt.subplots(ncols=2,nrows=1, figsize=(20, 15), sharey=True)

# axs.tick_params(axis='x', labelrotation=90)
# get unique variant types

var_con_scores_pre = sns.boxplot(y=exomiser_csv['score'],x=exomiser_csv['variant_consequence'],
                    color=".8", linewidth=.75, ax=axs[0]).set(title = "Exomiser scores across variant consequences",
                                                              ylabel = "Exomiser score",
                                                              xlabel = "Variant consequence")

var_con_scores_post = sns.boxplot(y=filtered_csv['score'],x=filtered_csv['variant_consequence'],
                     color=".8", linewidth=.75, ax=axs[1]).set(title = "Exomiser scores across variant consequences",
                                                              ylabel = "Exomiser score",
                                                              xlabel = "Variant consequence")
#var_con_scores_post.set_xticklabels(var_con_scores_post.get_xticklabels(),rotation=90)

# axs.set_xticks(ax.get_xticks(), axs.get_xticklabels(), rotation=45, ha='right')

fig.autofmt_xdate(rotation=90)
fig.savefig('all_variant_consequence_dev_3.png', bbox_inches ="tight", pad_inches = 0.5)


# violin plot of different exomiser scores of  different variant consequences 

fig, axs = plt.subplots(ncols=2,nrows=1, figsize=(20, 15), sharey=True)

var_clinvar_scores_pre = sns.boxplot(y=exomiser_csv['score'],x=exomiser_csv['clinvar'],
                      color=".8", linewidth=.75, ax=axs[0]).set(title = "Exomiser scores across clinvar classes",
                                                              ylabel = "Exomiser score",
                                                              xlabel = "Clinvar class")
#var_clinvar_scores_pre.set_xticklabels(var_clinvar_scores_pre.get_xticklabels(),rotation=90)

var_clinvar_scores_post = sns.boxplot(y=filtered_csv['score'],x=filtered_csv['clinvar'],
                      color=".8", linewidth=.75, ax=axs[1]).set(title = "Exomiser scores across clinvar classes",
                                                              ylabel = "Exomiser score",
                                                              xlabel = "Clinvar class")
#var_clinvar_scores_post.set_xticklabels(var_clinvar_scores_post.get_xticklabels(),rotation=90)

fig.autofmt_xdate(rotation=90)
fig.savefig('all_clinvar_annot_dev_3.png', bbox_inches ="tight", pad_inches = 0.5)
# fig.savefig("clinvar_post_filtering_output.png",  bbox_inches ="tight", pad_inches = 1)
