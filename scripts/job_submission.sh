#BSUB -q short
#BSUB -P gecip_lsf_access # This code will need to be changed to reflect your LSF project code. Please see the user guide for the list of registered codes.
#BSUB -o logs/%J.stdout
#BSUB -e logs/%J.stderr
#BSUB -cwd .
#BSUB -n 4
#BSUB -R "rusage[mem=36000] span[hosts=1]"
#BSUB -J Find_variants
export TMPDIR=/re_scratch/re_gecip/enhanced_interpretation/emackenzie

# set up environment

source /resources/conda/miniconda3/bin/activate
conda activate py3pypirev2

module purge
module load singularity/4.1.1
module load bcftools/1.16
module load samtools/1.16.1

mkdir -p ./logs

datetime=$(date +%y%m%d%H%M)


# run process over genes of interest

while read GENE_OF_INTEREST; do 

	outputfolder="${GENE_OF_INTEREST}""_"$datetime
	mkdir "$outputfolder"
	cd $outputfolder

	python3 ../query_exomiser.py "${GENE_OF_INTEREST}"

	bash ../run_vep.sh "${GENE_OF_INTEREST}"

	python3 ../merge_vcf_csv.py "${GENE_OF_INTEREST}"

	cd ..
done < genes_of_interest.txt
