#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH -A ACF-UTK0011
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --ntasks-per-node=12
#SBATCH --job-name=Decontam
#SBATCH --output=jobout-Decontam.o%j.txt
#SBATCH --chdir=/lustre/isaac24/proj/UTK0406/EmilyCantrell
#SBATCH --mail-user=sbn572@vols.utk.edu
#SBATCH --mail-type=ALL


# activate conda 
module load anaconda3/2024.06
source activate /lustre/isaac24/proj/UTK0406/shared_conda_envs/qiime2-amplicon-2025.7

THREADS=12

echo "Number of threads to use: $THREADS"


cd /lustre/isaac24/proj/UTK0406/EmilyCantrell

mkdir -p /lustre/isaac24/proj/UTK0406/EmilyCantrell/decontamination

#Decontamination scores - SR1
qiime quality-control decontam-identify \
    --i-table SR1_denoise/SR1_dada2_single_table.qza \
    --m-metadata-file /lustre/isaac24/proj/UTK0406/EmilyCantrell/NIJ_metadata.txt \
    --p-method prevalence \
    --p-prev-control-column treatment \
    --p-prev-control-indicator extraction_control \
    --o-decontam-scores decontamination/SR1_decontam_scores.qza

# SR1  score visualization
qiime quality-control decontam-score-viz \
    --i-decontam-scores decontamination/SR1_decontam_scores.qza \
    --i-table SR1_denoise/SR1_dada2_single_table.qza \
    --i-rep-seqs SR1_denoise/SR1_dada2_single_repseqs.qza \
    --p-threshold 0.1 \
    --p-no-weighted \
    --p-bin-size 0.05 \
    --o-visualization decontamination/SR1_decontam_scores.qzv


# Remove contaminants - SR1
qiime feature-table filter-features \
    --i-table SR1_denoise/SR1_dada2_single_table.qza \
    --m-metadata-file decontamination/SR1_decontam_scores.qza \
    --p-where '[p] >0.1 OR [p] IS NULL' \
    --o-filtered-table decontamination/SR1_filter_decontam_table.qza


# visualize rep seqs
qiime feature-table filter-seqs \
    --i-data SR1_denoise/SR1_dada2_single_repseqs.qza \
    --i-table decontamination/SR1_filter_decontam_table.qza \
    --o-filtered-data decontamination/SR1_filter_rep_seq.qza

# visualize table
qiime feature-table summarize \
    --i-table decontamination/SR1_filter_decontam_table.qza \
    --o-visualization decontamination/SR1_filter_decontam_table.qzv




# Decontamination scores - SR2
qiime quality-control decontam-identify \
    --i-table SR2_denoise/SR2_dada2_single_table.qza \
    --m-metadata-file /lustre/isaac24/proj/UTK0406/EmilyCantrell/NIJ_metadata.txt \
    --p-method prevalence \
    --p-prev-control-column treatment \
    --p-prev-control-indicator extraction_control \
    --o-decontam-scores decontamination/SR2_decontam_scores.qza 

# SR2 decontam score visualization
qiime quality-control decontam-score-viz \
    --i-decontam-scores decontamination/SR2_decontam_scores.qza \
    --i-table SR2_denoise/SR2_dada2_single_table.qza \
    --i-rep-seqs SR2_denoise/SR2_dada2_single_repseqs.qza \
    --p-threshold 0.1 \
    --p-no-weighted \
    --p-bin-size 0.05 \
    --o-visualization decontamination/SR2_decontam_scores.qzv


# Remove contaminants - SR2
qiime feature-table filter-features \
    --i-table SR2_denoise/SR2_dada2_single_table.qza \
    --m-metadata-file decontamination/SR2_decontam_scores.qza \
    --p-where '[p] >0.1 OR [p] IS NULL' \
    --o-filtered-table decontamination/SR2_filter_decontam_table.qza


# visualize rep seqs
qiime feature-table filter-seqs \
    --i-data SR2_denoise/SR2_dada2_single_repseqs.qza \
    --i-table decontamination/SR2_filter_decontam_table.qza \
    --o-filtered-data decontamination/SR2_filter_rep_seq.qza


#visualize table
qiime feature-table summarize \
    --i-table decontamination/SR2_filter_decontam_table.qza \
    --o-visualization decontamination/SR2_filter_decontam_table.qzv