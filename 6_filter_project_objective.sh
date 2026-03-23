#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -A ACF-UTK0011
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --ntasks-per-node=12
#SBATCH --job-name=filter
#SBATCH --output=jobout-filter.o%j.txt
#SBATCH --chdir=/lustre/isaac24/proj/UTK0406/EmilyCantrell
#SBATCH --mail-user=sbn572@vols.utk.edu
#SBATCH --mail-type=ALL


# activate conda 
module load anaconda3/2024.06
source activate /lustre/isaac24/proj/UTK0406/shared_conda_envs/qiime2-amplicon-2025.4

THREADS=12

echo "Number of threads to use: $THREADS"

cd /lustre/isaac24/proj/UTK0406/EmilyCantrell


# filter data by proj objective

# Trace Evidence Soil
qiime feature-table filter-samples \
--i-table taxonomy/SR1_SR2_tax_filter.qza \
--m-metadata-file NIJ_metadata.txt \
--p-where "[project_objective]='Trace_Evidence_Soil'" \
--o-filtered-table merged/TE_soil_table.qza


# visualization
qiime feature-table summarize \
    --i-table merged/TE_soil_table.qza \
    --m-sample-metadata-file NIJ_metadata.txt \
    --o-visualization merged/TE_soil_table.qzv



# Trace Evidence Touch
qiime feature-table filter-samples \
--i-table taxonomy/SR1_SR2_tax_filter.qza \
--m-metadata-file NIJ_metadata.txt \
--p-where "[project_objective]='Trace_Evidence_Touch'" \
--o-filtered-table merged/TE_touch_table.qza


# visualization
qiime feature-table summarize \
    --i-table merged/TE_touch_table.qza \
    --m-sample-metadata-file NIJ_metadata.txt \
    --o-visualization merged/TE_touch_table.qzv



# Trace evidence touch - pvc
qiime feature-table filter-samples \
    --i-table merged/TE_touch_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location] = 'left_hand_pvc' OR [sample_location] = 'right_hand_pvc'" \
    --o-filtered-table merged/TE_pvc_table.qza

# visualization
qiime feature-table summarize \
    --i-table merged/TE_pvc_table.qza \
    --m-sample-metadata-file NIJ_metadata.txt \
    --o-visualization merged/TE_pvc_table.qzv


# Trace evidence touch - hand
qiime feature-table filter-samples \
    --i-table merged/TE_touch_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location] = 'left_hand_palm' OR [sample_location] = 'right_hand_palm'" \
    --o-filtered-table merged/TE_palm_table.qza

# visualization
qiime feature-table summarize \
    --i-table merged/TE_palm_table.qza \
    --m-sample-metadata-file NIJ_metadata.txt \
    --o-visualization merged/TE_palm_table.qzv
