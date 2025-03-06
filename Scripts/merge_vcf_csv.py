
import sys
import numpy as np
import functools
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import vcf

# get gene name
print("Merging CSV and VCF outputs for gene: ", sys.argv[1])
gene_name = sys.argv[1]

# read in vcf and csv outputs from query_exomiser.py and run_vep.sh
exomiser_output_name = gene_name + "_output_table.csv"
exomiser_csv = pd.read_csv(exomiser_output_name)
exomiser_csv = exomiser_csv.drop("Unnamed: 0", axis = 1)
exomiser_csv =  exomiser_csv.dropna(subset=['hgvs'])
vcf_output_name = gene_name + "_vep_input_annotated.vcf"

# create empty columns for adding in annotations
exomiser_csv = exomiser_csv.assign(variant_type="NA", variant_consequence="NA",
	existing_variant="NA", clinvar="NA", clinvar_evidence="NA",
	gnomad=" ", REVEL=" ", CADD_PHRED=" ")

# print(exomiser_csv)

vcf_reader = vcf.Reader(open(vcf_output_name, 'r'))
for record in vcf_reader:
	# this is a dictionary
	info = record.INFO
	CSQ = info["CSQ"]
	# split the CSQ entry into individual entries, which have constant positions
        # from CSQ want VARIANT_CLASS (3), consequnece (4), feature (8), existing variant (12),
        # clinvar clin sig (15, clinvar cin sig conf (16), gnomad (18), REVEL(57), CADD_PHRED (58)
	VEP_info_Split = str(CSQ[0]).split("|")

	# find corresponding line 
	HGVS = list(info.keys())[0]
	# inconsistent format in p. of hgvs for different variants, so split
	if "*" in HGVS:
		HGVS_split = HGVS.split("p.")
		HGVS = HGVS_split[0]

		index = exomiser_csv[exomiser_csv['hgvs'].str.contains(HGVS, regex=False)].index.tolist()
		exomiser_csv.at[index, "variant_type"]=VEP_info_Split[3]
		exomiser_csv.at[index, "variant_consequence"]=VEP_info_Split[4]
		exomiser_csv.at[index, "existing_variant"]=VEP_info_Split[12]
		exomiser_csv.at[index, "clinvar"]=VEP_info_Split[15]
		exomiser_csv.at[index, "clinvar_evidence"]=VEP_info_Split[16]
		exomiser_csv.at[index, "gnomad"]=VEP_info_Split[18]
		exomiser_csv.at[index, "REVEL"]=VEP_info_Split[57]
		exomiser_csv.at[index, "CADD_PHRED"]=VEP_info_Split[58]

	elif "+" in HGVS:

		index = exomiser_csv[exomiser_csv['hgvs'] == HGVS].index.tolist()
		exomiser_csv.at[index, "variant_type"]=VEP_info_Split[3]
		exomiser_csv.at[index, "variant_consequence"]=VEP_info_Split[4]
		exomiser_csv.at[index, "existing_variant"]=VEP_info_Split[12]
		exomiser_csv.at[index, "clinvar"]=VEP_info_Split[15]
		exomiser_csv.at[index, "clinvar_evidence"]=VEP_info_Split[16]
		exomiser_csv.at[index, "gnomad"]=VEP_info_Split[18]
		exomiser_csv.at[index, "REVEL"]=VEP_info_Split[57]
		exomiser_csv.at[index, "CADD_PHRED"]=VEP_info_Split[58]
	else:
		HGVS_split = HGVS.split("p.")
		HGVS = HGVS_split[0]

		index = exomiser_csv[exomiser_csv['hgvs'].str.contains(HGVS)].index.tolist()
		exomiser_csv.at[index, "variant_type"]=VEP_info_Split[3]
		exomiser_csv.at[index, "variant_consequence"]=VEP_info_Split[4]
		exomiser_csv.at[index, "existing_variant"]=VEP_info_Split[12]
		exomiser_csv.at[index, "clinvar"]=VEP_info_Split[15]
		exomiser_csv.at[index, "clinvar_evidence"]=VEP_info_Split[16]
		exomiser_csv.at[index, "gnomad"]=VEP_info_Split[18]
		exomiser_csv.at[index, "REVEL"]=VEP_info_Split[57]
		exomiser_csv.at[index, "CADD_PHRED"]=VEP_info_Split[58]


## filtering stage ##


# remove those with computational evidence prediciting benign
# remove clinvar bengin variants
# remove variants in non-coding regions

# print(exomiser_csv)
csv_name = gene_name + "_merged_output.csv"
filtered_name = gene_name + "_filtered_VUS.csv"
participants = gene_name + "_participants_in_filtered.csv"

exomiser_csv.to_csv(csv_name, index=False)

# remove benign 
exomiser_csv = exomiser_csv[~exomiser_csv['clinvar'].str.contains('benign')]
exomiser_csv = exomiser_csv[~exomiser_csv['variant_consequence'].isin(['synonymous_variant', 'intron_variant',
               'upstream_gene_variant', 'downstream_gene_variant', 'intergenic_variant', '5_prime_UTR_variant',
               '3_prime_UTR_variant'])]

exomiser_csv = exomiser_csv[~exomiser_csv['clinvar'].str.contains('Benign')]

# remove low CADD, large indels have no CADD
exomiser_csv['CADD_PHRED'] = pd.to_numeric(exomiser_csv['CADD_PHRED'], errors='coerce')
exomiser_csv = exomiser_csv[ (exomiser_csv['CADD_PHRED'] >= 22) | (exomiser_csv['CADD_PHRED'].isnull()) ]


# remove high sift
exomiser_csv['sift'] = pd.to_numeric(exomiser_csv['sift'], errors='coerce')
exomiser_csv  = exomiser_csv[ (exomiser_csv['sift'].isnull()) | (exomiser_csv['sift'] < 0.327) ]

# remove low polyphen
exomiser_csv['poly_phen'] = pd.to_numeric(exomiser_csv['poly_phen'], errors='coerce')
exomiser_csv  = exomiser_csv[ (exomiser_csv['poly_phen'].isnull()) | (exomiser_csv['poly_phen'] > 0.113) ]

# remove low REVEL
exomiser_csv['REVEL'] = pd.to_numeric(exomiser_csv['REVEL'], errors='coerce')
exomiser_csv  = exomiser_csv[ (exomiser_csv['REVEL'] > 0.5) | (exomiser_csv['REVEL'].isnull()) ]


# Get participant IDs for airlock
partic_ID = exomiser_csv['participant_id'].to_list()
# remove paricipant IDs from table
exomiser_csv_no_part = exomiser_csv.drop('participant_id', axis=1)

# create filtered csv
exomiser_csv_no_part.to_csv(filtered_name, index=False)

# create file of participants in filtered csv
with open(participants, 'w') as file:
	for item in partic_ID:
	# write each item on a new line
		file.write("%s\n" % item)



# make figures of pre and post filtering

# violin plot of exomiser score different variant types

# set up dataframe 

#fig, axs = plt.subplots(ncols=2,nrows=1, figsize=(8, 10))

# get unique variant types

#var_con_scores_pre = sns.violinplot(y=exomiser_csv['score'],x=exomiser_csv['variant_consequence'], widths=0.1,
#                     showmeans=True, showmedians=True, inner="point", cut=0, ax=axs[0])
#var_con_scores_pre.tick_params(axis='x', labelrotation=90)#
#
#var_con_scores_post = sns.violinplot(y=high_revel['score'],x=high_revel['variant_consequence'], widths=0.1, 
#                     showmeans=True, showmedians=True, inner="point", cut=0, ax=axs[1])
#var_con_scores_post.tick_params(axis='x', labelrotation=90)
#
#fig.savefig('variant_consequence.png', bbox_inches ="tight", pad_inches = 1)


# violin plot of different exomiser scores of  different variant consequences 

#fig, axs = plt.subplots(ncols=2,nrows=1, figsize=(4, 4))

#var_clinvar_scores_pre = sns.violinplot(y=exomiser_csv['score'],x=exomiser_csv['clinvar'],
#                     showmeans=True, showmedians=True, inner="point", cut=0, ax=axs[0])
#var_clinvar_scores_pre.tick_params(axis='x', labelrotation=90)
#
#var_clinvar_scores_post = sns.violinplot(y=high_revel['score'],x=high_revel['clinvar'],
#                     showmeans=True, showmedians=True, inner="point", cut=0, ax=axs[1])
#var_clinvar_scores_post.tick_params(axis='x', labelrotation=90)

#fig.savefig('clinvar_annot.png', bbox_inches ="tight", pad_inches = 1)
# fig.savefig("clinvar_post_filtering_output.png",  bbox_inches ="tight", pad_inches = 1)

