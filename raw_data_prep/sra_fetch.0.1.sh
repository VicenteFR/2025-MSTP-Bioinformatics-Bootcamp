# ---> Retrieve public data from SRA's bioproject PRJNA335583
# Bioproject from paper titled: Drug-tolerant persister cancer cells are vulnerable to GPX4 inhibition
# RunInfo and Accession List retrieved from this page: https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=335583

# ---> Arguments.
# General files.
BIOPRJ=PRJNA335583
# Path and file definition.
SAMPLE_LIST_FILE=/tscc/nfs/home/${USER}/raw_data/SraAccList.csv
R_SCRIPT=/tscc/nfs/home/${USER}/raw_data/jobs_scripts/sra_rename.0.1.R
GEN_OUT_DIR=/tscc/lustre/ddn/scratch/${USER}/raw_data/${BIOPRJ}
if [ ! -d ${GEN_OUT_DIR} ]; then mkdir -p ${GEN_OUT_DIR}; fi
OUT_DIR=${GEN_OUT_DIR}/fastqs
if [ ! -d ${OUT_DIR} ]; then mkdir ${OUT_DIR}; fi
JOB_DIR=${GEN_OUT_DIR}/jobs_scripts
if [ ! -d ${JOB_DIR} ]; then mkdir ${JOB_DIR}; fi

# ---> Download dataset.
# Fetch all data.
# cd ${OUT_DIR}
# prefetch --option-file ${SAMPLE_LIST_FILE}
# Dump fastqs.
# for SAMPLE in $( cat ${SAMPLE_LIST_FILE} ); do
#     TMP_PATH=${OUT_DIR}/${SAMPLE}
#     cd ${TMP_PATH}
#     # pwd
#     fasterq-dump --split-files ${SAMPLE}.sra
# done
# Fetch as fastq files in a single step.
for SAMPLE in $( cat ${SAMPLE_LIST_FILE} ); do
    # Define files.
    JOB_FILE=${JOB_DIR}/${SAMPLE}.job.sh
    OUT_FILE=${JOB_DIR}/${SAMPLE}.out.txt
    ERR_FILE=${JOB_DIR}/${SAMPLE}.err.txt
    cat <<-EOF > ${JOB_FILE}
#!/bin/bash
#SBATCH --job-name=${SAMPLE}_fetch
#SBATCH --output=${OUT_FILE}
#SBATCH --error=${ERR_FILE}
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=4G
#SBATCH --time=3:00:00
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --account=htl191
cd ${OUT_DIR}
fasterq-dump --split-files ${SAMPLE}
EOF
    sbatch ${JOB_FILE}
done

# ---> Rename files on the basis of their SRR IDs.
JOB_FILE=${JOB_DIR}/rename_files.job.sh
OUT_FILE=${JOB_DIR}/rename_files.out.txt
ERR_FILE=${JOB_DIR}/rename_files.err.txt
cat <<-EOF > ${JOB_FILE}
#!/bin/bash
#SBATCH --job-name=rename_files
#SBATCH --output=${OUT_FILE}
#SBATCH --error=${ERR_FILE}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --account=htl191
conda init bash
source ~/.bashrc
conda activate 2025-mstp-bootcamp
Rscript ${R_SCRIPT}
EOF
sbatch ${JOB_FILE}

# ---> Fastq basic stats
JOB_FILE=${JOB_DIR}/basic_stats.job.sh
OUT_FILE=${JOB_DIR}/basic_stats.out.txt
ERR_FILE=${JOB_DIR}/basic_stats.err.txt
cat <<-EOF > ${JOB_FILE}
#!/bin/bash
#SBATCH --job-name=basic_stats
#SBATCH --output=${OUT_FILE}
#SBATCH --error=${ERR_FILE}
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --account=htl191
seqkit stats \
--threads 2 \
--tabular ${OUT_DIR}/*.fastq \
> ${OUT_DIR}/fastq_stats.tsv
EOF
sbatch ${JOB_FILE}

# ---> Zip fastq files.
JOB_FILE=${JOB_DIR}/zip_files.job.sh
OUT_FILE=${JOB_DIR}/zip_files.out.txt
ERR_FILE=${JOB_DIR}/zip_files.err.txt
cat <<-EOF > ${JOB_FILE}
#!/bin/bash
#SBATCH --job-name=zip_files
#SBATCH --output=${OUT_FILE}
#SBATCH --error=${ERR_FILE}
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=4G
#SBATCH --time=5:00:00
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --account=htl191
gzip ${OUT_DIR}/*.fastq
EOF
sbatch ${JOB_FILE}