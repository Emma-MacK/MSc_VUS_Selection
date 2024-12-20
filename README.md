# MSc_VUS_Selection #
A repository holding code used to query tables in the GEL research environment and select VUS for functional studies
This code was developed as part of an MSc project aiming to identify VUS that would benefit from a functional study. 
The code was designed to run within the GEL research environment.
To run this code you must have acces to and be within the GEL research environment

## What are the scripts ##

### job_submission.sh ###

This script is used to run query_exomiser.py, run_vep.sh, and merge_vcf_csv.py as a submitted job on the GEL research environment HPC
Requires:
* A txt file containing the symbols of genes to be processed, with each gene symbol on a seperate line

### query_exomiser.py ###

This script runs an SQL command via the python package LabKey to select variants in a gene of interest from the GEL research environment exomiser table. Variants are only selected if the participant ID corresponds to a participant that is a proband, has consented for research, and has not had a submission to the diagnostic_discovery table. The names of panels a participant was tested for is extracted from the panels_applied table. Chromosome, position, reference, alternate, and HGVS are extracted to build a VCF file to be used by the run_vep.sh script. The script outputs a table of variant information, having one unique entry for each HGVS ID. The row that had the highest exomiser score is kept when condensing the table down to one entry per HGVS. Output files are saved with the gene symbol in the file name.

Requires:
* An gene symbol
 
Outputs:
* *_vep_input.txt - a file in VCF format to pass variants to run_vep.sh
* *_output_table.csv - a table of variants showing the entry with the top exomiser score for each HGVS

### run_vep.sh ###

This script creates links to VEP annotation sources, and runs VEP on the given *_vep_input.txt file produced by query_exomiser.py. This outputs an annotated vcf file.

Requires:
* *_vep_input.txt

Outputs:
* *_annotated.vcf

### merge_vcf_csv.py ###

This script combines the *_output_table.csv produced by query_exomiser.py and the *_annotated.vcf file produced by run_vep.sh to keep specific annotations of interest (clinvar annotation, consequence, CADD_PHRED, REVEL). A csv file is produced of the combined table but would not be able to be extracted from the airlock. This script then performs filtering, removing variants that fall into non-coding regions or where computational evidence is present indicating a variant is benign. The participants ID are removed from the filtered table and added to their own file. This is done as Participant IDs cannot be extracted, but the airlock process requires the IDs of Participants who's data is being exported.

Requires:
* *_annotated.vcf
* *_output_table.csv

Outputs:
* *_merged_output.csv
* *_filtered_VUS.csv
* *_participants_in_filtered.csv

### create_figures.py ###

This script is used to analyse the range of exomiser scores across both consequence type and clinvar predictions pre and post filtering. It is run seperately from the other scripts and requires a file containing all of the variants across given genes of interest before filtering and a file containing all of the variants across genes of interest after filtering

Requires:
* all_genes_merged_variants.csv
* all_genes_filtered_variants_stage_4.csv - currently this file has "stage_4" in the name due to it being the fourth "stage" of developing this pipeline

Outputs:
* all_variant_consequence_dev_3.png
* all_clinvar_annot_dev_3.png

## Set up the GEL research environment ##



## Submitting to the HPC ##



## Running scripts individually ## 
