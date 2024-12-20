
GENE_OF_INTEREST=$1

# Set variables
IMG=/gel_data_resources/containers/singularity_containers/VEP/v109/vep_109.sif
CACHE=/public_data_resources/vep_resources/VEP_109/
REFFASTA=/public_data_resources/reference/GRCh38/GRCh38Decoy_no_alt.fa
INPUT_FILE="${GENE_OF_INTEREST}"_vep_input.txt # This file will need to be created within your working directory

sort -k 2 "${INPUT_FILE}" > input.txt

# Bind paths
# some paths are specific to my work environment
MOUNT_GENOMES='/nas/weka.gel.zone/pgen_genomes:/nas/weka.gel.zone/pgen_genomes,/genomes:/genomes'
MOUNT_GEL_DATA_RESOURCES='/nas/weka.gel.zone/pgen_int_data_resources:/nas/weka.gel.zone/pgen_int_data_resources,/gel_data_resources:/gel_data_resources'
MOUNT_PUBLIC_DATA_RESOURCES='/nas/weka.gel.zone/pgen_public_data_resources:/nas/weka.gel.zone/pgen_public_data_resources,/public_data_resources:/public_data_resources'
MOUNT_SCRATCH='/nas/weka.gel.zone/re_scratch:/nas/weka.gel.zone/re_scratch,/re_scratch:/re_scratch'
MOUNT_WD='/nas/weka.gel.zone/home/emackenzie/VEP_test/109:/nas/weka.gel.zone/home/emackenzie/VEP_test/109'

# Set plugin variables
CADD16='/public_data_resources/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz'
REVEL='/public_data_resources/vep_resources/REVEL/revel_v1.3_GRCh38.tsv.gz'
GNOMAD='/public_data_resources/vep_resources/VEP_plugins/gnomADc.pm'
SPLICEAIRAW38='/public_data_resources/SpliceAI/Predicting_splicing_from_primary_sequence-66029966/genome_scores_v1.3/spliceai_scores.raw.snv.hg38.vcf.gz'
SPLICEAIINDEL38='/public_data_resources/SpliceAI/Predicting_splicing_from_primary_sequence-66029966/genome_scores_v1.3/spliceai_scores.raw.indel.hg38.vcf.gz'

mkdir -p ./logs

output="$(basename "${INPUT_FILE}" | cut -f 1 -d '.')"
echo "${output}"

# Run VEP 109 
singularity exec --bind "${MOUNT_WD}","${MOUNT_GENOMES}","${MOUNT_GEL_DATA_RESOURCES}","${MOUNT_PUBLIC_DATA_RESOURCES}","${MOUNT_SCRATCH}" "${IMG}" vep \
    --cache \
    --species homo_sapiens \
    --vcf \
    --pick \
    --force_overwrite \
    --assembly GRCh38 \
    --symbol --hgvs --hgvsg --check_existing --variant_class \
    --plugin CADD,"${CADD16}" \
    --plugin SpliceAI,snv="${SPLICEAIRAW38}",indel="${SPLICEAIINDEL38}" \
    --plugin REVEL,"${REVEL}" \
    --plugin gnomADc,"${GNOMAD}" \
    --dir_cache "${CACHE}" \
    --fasta "${REFFASTA}" \
    --input_file input.txt \
    --output_file "${PWD}"/"${output}"_annotated.vcf \
    --no_check_variants_order \
    --custom /public_data_resources/clinvar/20240627/vcf_GRCh38/clinvar_20240624.vcf.gz,ClinVar,vcf,exact,0,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNVC,CLNVCSO,CLNVI \
    --offline \
    --fields Allele,SYMBOL,HGNC_ID,VARIANT_CLASS,Consequence,IMPACT,EXON,INTRON,Feature,HGVSc,HGVSp,HGVS_OFFSET,Existing_variation,STRAND,ClinVar,ClinVar_CLNSIG,ClinVar_CLNSIGCONF,ClinVar_CLNDN,gnomADc_AC,gnomADg_AN,gnomADg_AF,gnomADg_nhomalt,gnomADg_popmax,gnomADg_AC_popmax,gnomADg_AN_popmax,gnomADg_AF_popmax,gnomADg_nhomalt_popmax,gnomADe_AC,gnomADe_AN,gnomADe_AF,gnomADe_nhomalt,gnomADe_popmax,gnomADe_AC_popmax,gnomADe_AN_popmax,gnomADe_AF_popmax,gnomADe_nhomalt_popmax,gnomADe_non_cancer_AC,gnomADe_non_cancer_AN,gnomADe_non_cancer_AF,gnomADe_non_cancer_nhomalt,gnomADe_non_cancer_AC_popmax,gnomADe_non_cancer_AN_popmax,gnomADe_non_cancer_AF_popmax,gnomADe_non_cancer_nhomalt_popmax,gnomADe_non_cancer_popmax,HGMD,HGMD_PHEN,HGMD_CLASS,HGMD_RANKSCORE,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL,SpliceAI_pred_DP_AG,SpliceAI_pred_DP_AL,SpliceAI_pred_DP_DG,SpliceAI_pred_DP_DL,REVEL,CADD_PHRED 


