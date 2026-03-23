#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -A ACF-UTK0011
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --ntasks-per-node=12
#SBATCH --job-name=rarefaction
#SBATCH --output=jobout-rarefaction.o%j.txt
#SBATCH --chdir=/lustre/isaac24/proj/UTK0406/EmilyCantrell
#SBATCH --mail-user=sbn572@vols.utk.edu
#SBATCH --mail-type=ALL


# activate conda 
module load anaconda3/2024.06
source activate /lustre/isaac24/proj/UTK0406/shared_conda_envs/qiime2-amplicon-2025.7

THREADS=12

echo "Number of threads to use: $THREADS"

cd /lustre/isaac24/proj/UTK0406/EmilyCantrell

mkdir -p rarefaction

# all samples
qiime diversity alpha-rarefaction \
--i-table taxonomy/SR1_SR2_tax_filter.qza \
--i-phylogeny tree/tree.qza \
--m-metadata-file NIJ_metadata.txt \
--p-min-depth 10 \
--p-max-depth 10000 \
--o-visualization rarefaction/alpha_rarefaction_curves.qzv \


# TE - soil
qiime diversity alpha-rarefaction \
    --i-table merged/TE_soil_table.qza \
    --i-phylogeny tree/tree.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-min-depth 10 \
    --p-max-depth 25000 \
    --o-visualization rarefaction/TE_soil_rarefaction.qzv


# TE - skin
qiime diversity alpha-rarefaction \
    --i-table merged/TE_touch_table.qza \
    --i-phylogeny tree/tree.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-min-depth 10 \
    --p-max-depth 25000 \
    --o-visualization rarefaction/TE_touch_rarefaction.qzv




