#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -A ACF-UTK0011
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --ntasks-per-node=12
#SBATCH --job-name=Merge
#SBATCH --output=jobout-Merge.o%j.txt
#SBATCH --chdir=/lustre/isaac24/proj/UTK0406/EmilyCantrell
#SBATCH --mail-user=sbn572@vols.utk.edu
#SBATCH --mail-type=ALL


# activate conda 
module load anaconda3/2024.06
source activate /lustre/isaac24/proj/UTK0406/shared_conda_envs/qiime2-amplicon-2025.7

THREADS=12

echo "Number of threads to use: $THREADS"


cd /lustre/isaac24/proj/UTK0406/EmilyCantrell


mkdir -p /lustre/isaac24/proj/UTK0406/EmilyCantrell/merged


# merge tables
qiime feature-table merge \
    --i-tables decontamination/SR1_filter_decontam_table.qza \
    --i-tables decontamination/SR2_filter_decontam_table.qza \
    --o-merged-table merged/SR1_SR2_merged_table.qza


# merge seqs
qiime feature-table merge-seqs \
    --i-data decontamination/SR1_filter_rep_seq.qza \
    --i-data decontamination/SR2_filter_rep_seq.qza \
    --o-merged-data merged/SR1_SR2_merged_rep_seqs.qza


# visualize table
qiime feature-table summarize \
    --i-table merged/SR1_SR2_merged_table.qza \
    --o-visualization merged/SR1_SR2_merged_table.qzv \
    --m-sample-metadata-file NIJ_metadata.txt


# visualize seqs
qiime feature-table tabulate-seqs \
    --i-data merged/SR1_SR2_merged_rep_seqs.qza \
    --o-visualization merged/SR1_SR2_merged_rep_seqs.qzv