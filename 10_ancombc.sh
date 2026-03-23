#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -A ACF-UTK0011
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --ntasks-per-node=12
#SBATCH --job-name=ancombc
#SBATCH --output=jobout-ancombc.o%j.txt
#SBATCH --chdir=/lustre/isaac24/proj/UTK0406/EmilyCantrell
#SBATCH --mail-user=sbn572@vols.utk.edu
#SBATCH --mail-type=ALL


# activate conda 
module load anaconda3/2024.06
source activate /lustre/isaac24/proj/UTK0406/shared_conda_envs/qiime2-amplicon-2025.4
THREADS=12

echo "Number of threads to use: $THREADS"

cd /lustre/isaac24/proj/UTK0406/EmilyCantrell



### Filter table by location

echo "filter locations"

# Filter location A
 qiime feature-table filter-samples \
    --i-table merged/TE_soil_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location] = 'location_A'" \
    --o-filtered-table new-TE-soil-coremetrics/locationA_table.qza



# Filter location B
 qiime feature-table filter-samples \
    --i-table merged/TE_soil_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location] = 'location_B'" \
    --o-filtered-table new-TE-soil-coremetrics/locationB_table.qza



# Filter location C
 qiime feature-table filter-samples \
    --i-table merged/TE_soil_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location] = 'location_C'" \
    --o-filtered-table new-TE-soil-coremetrics/locationC_table.qza



# Filter location D
 qiime feature-table filter-samples \
    --i-table merged/TE_soil_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location] = 'location_D'" \
    --o-filtered-table new-TE-soil-coremetrics/locationD_table.qza




echo "start ancombc2"

# location A ancombc2
qiime composition ancombc2 \
    --i-table new-TE-soil-coremetrics/locationA_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-fixed-effects-formula treatment \
    --p-reference-levels treatment::control \
    --o-ancombc2-output new-TE-soil-coremetrics/locationA_ancombc2.qza 


# location B ancombc2
qiime composition ancombc2 \
    --i-table new-TE-soil-coremetrics/locationB_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-fixed-effects-formula treatment \
    --p-reference-levels treatment::control \
    --o-ancombc2-output new-TE-soil-coremetrics/locationB_ancombc2.qza 


# location C ancombc2
qiime composition ancombc2 \
    --i-table new-TE-soil-coremetrics/locationC_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-fixed-effects-formula treatment \
    --p-reference-levels treatment::control \
    --o-ancombc2-output new-TE-soil-coremetrics/locationC_ancombc2.qza 


# location D ancombc2
qiime composition ancombc2 \
    --i-table new-TE-soil-coremetrics/locationD_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-fixed-effects-formula treatment \
    --p-reference-levels treatment::control \
    --o-ancombc2-output new-TE-soil-coremetrics/locationD_ancombc2.qza 



echo "start visualizations"

# location A abc2
qiime composition ancombc2-visualizer \
    --i-data new-TE-soil-coremetrics/locationA_ancombc2.qza \
    --i-taxonomy taxonomy/taxonomy.qza \
    --o-visualization new-TE-soil-coremetrics/locationA_ancombc2.qzv


# location B abc2
qiime composition ancombc2-visualizer \
    --i-data new-TE-soil-coremetrics/locationB_ancombc2.qza \
    --i-taxonomy taxonomy/taxonomy.qza \
    --o-visualization new-TE-soil-coremetrics/locationB_ancombc2.qzv


# location C abc2
qiime composition ancombc2-visualizer \
    --i-data new-TE-soil-coremetrics/locationC_ancombc2.qza \
    --i-taxonomy taxonomy/taxonomy.qza \
    --o-visualization new-TE-soil-coremetrics/locationC_ancombc2.qzv


# location D abc2
qiime composition ancombc2-visualizer \
    --i-data new-TE-soil-coremetrics/locationD_ancombc2.qza \
    --i-taxonomy taxonomy/taxonomy.qza \
    --o-visualization new-TE-soil-coremetrics/locationD_ancombc2.qzv


# export taxonomy
qiime tools export \
    --input-path taxonomy/taxonomy.qza \
    --output-path taxonomy/taxonomy_export


## export location files
qiime tools export \
    --input-path new-TE-soil-coremetrics/locationA_ancombc2.qza \
    --output-path new-TE-soil-coremetrics/locationA_export


qiime tools export \
    --input-path new-TE-soil-coremetrics/locationB_ancombc2.qza \
    --output-path new-TE-soil-coremetrics/locationB_export

qiime tools export \
    --input-path new-TE-soil-coremetrics/locationC_ancombc2.qza \
    --output-path new-TE-soil-coremetrics/locationC_export

qiime tools export \
    --input-path new-TE-soil-coremetrics/locationD_ancombc2.qza \
    --output-path new-TE-soil-coremetrics/locationD_export