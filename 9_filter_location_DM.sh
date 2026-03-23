#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -A ACF-UTK0011
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --ntasks-per-node=12
#SBATCH --job-name=filterDM
#SBATCH --output=jobout-filterDM.o%j.txt
#SBATCH --chdir=/lustre/isaac24/proj/UTK0406/EmilyCantrell
#SBATCH --mail-user=sbn572@vols.utk.edu
#SBATCH --mail-type=ALL


# activate conda 
module load anaconda3/2024.06
source activate /lustre/isaac24/proj/UTK0406/shared_conda_envs/qiime2-amplicon-2025.4
THREADS=12

echo "Number of threads to use: $THREADS"

cd /lustre/isaac24/proj/UTK0406/EmilyCantrell


echo "filter unweighted DM"

# filter unweighted distance  matrix - location A
qiime diversity filter-distance-matrix \
    --i-distance-matrix new-TE-soil-coremetrics/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_A'" \
    --o-filtered-distance-matrix new-TE-soil-coremetrics/locationA_unweighted_distance_matrix.qza


# filter unweighted distance  matrix - location B
qiime diversity filter-distance-matrix \
    --i-distance-matrix new-TE-soil-coremetrics/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_B'" \
    --o-filtered-distance-matrix new-TE-soil-coremetrics/locationB_unweighted_distance_matrix.qza

# filter unweighted distance matrix - location C
qiime diversity filter-distance-matrix \
    --i-distance-matrix new-TE-soil-coremetrics/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_C'" \
    --o-filtered-distance-matrix new-TE-soil-coremetrics/locationC_unweighted_distance_matrix.qza

 # filter unweighted distance matrix - location D
qiime diversity filter-distance-matrix \
    --i-distance-matrix new-TE-soil-coremetrics/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_D'" \
    --o-filtered-distance-matrix new-TE-soil-coremetrics/locationD_unweighted_distance_matrix.qza


echo "filter weighted DM"

# filter weighted DM location A
qiime diversity filter-distance-matrix \
    --i-distance-matrix new-TE-soil-coremetrics/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_A'" \
    --o-filtered-distance-matrix new-TE-soil-coremetrics/locationA_weighted_distance_matrix.qza


# filter weighted DM location B
qiime diversity filter-distance-matrix \
    --i-distance-matrix new-TE-soil-coremetrics/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_B'" \
    --o-filtered-distance-matrix new-TE-soil-coremetrics/locationB_weighted_distance_matrix.qza


# filter weighted DM location C
qiime diversity filter-distance-matrix \
    --i-distance-matrix new-TE-soil-coremetrics/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_C'" \
    --o-filtered-distance-matrix new-TE-soil-coremetrics/locationC_weighted_distance_matrix.qza


# filter weighted DM location D
qiime diversity filter-distance-matrix \
    --i-distance-matrix new-TE-soil-coremetrics/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_D'" \
    --o-filtered-distance-matrix new-TE-soil-coremetrics/locationD_weighted_distance_matrix.qza



echo " start unweighted pcoa"

# location A unweighted pcoa
qiime diversity pcoa \
    --i-distance-matrix new-TE-soil-coremetrics/locationA_unweighted_distance_matrix.qza \
    --o-pcoa new-TE-soil-coremetrics/locationA_unweighted_pcoa.qza


# location B unweighted pcoa
qiime diversity pcoa \
    --i-distance-matrix new-TE-soil-coremetrics/locationB_unweighted_distance_matrix.qza \
    --o-pcoa new-TE-soil-coremetrics/locationB_unweighted_pcoa.qza

# location C unweighted pcoa
qiime diversity pcoa \
    --i-distance-matrix new-TE-soil-coremetrics/locationC_unweighted_distance_matrix.qza \
    --o-pcoa new-TE-soil-coremetrics/locationC_unweighted_pcoa.qza

# location D unweighted pcoa
qiime diversity pcoa \
    --i-distance-matrix new-TE-soil-coremetrics/locationD_unweighted_distance_matrix.qza \
    --o-pcoa new-TE-soil-coremetrics/locationD_unweighted_pcoa.qza


echo "start weighted pcoa"


# location A weighted pcoa
qiime diversity pcoa \
    --i-distance-matrix new-TE-soil-coremetrics/locationA_weighted_distance_matrix.qza \
    --o-pcoa new-TE-soil-coremetrics/locationA_weighted_pcoa.qza


# location B weighted pcoa
qiime diversity pcoa \
    --i-distance-matrix new-TE-soil-coremetrics/locationB_weighted_distance_matrix.qza \
    --o-pcoa new-TE-soil-coremetrics/locationB_weighted_pcoa.qza 


# location C weighted pcoa
qiime diversity pcoa \
    --i-distance-matrix new-TE-soil-coremetrics/locationC_weighted_distance_matrix.qza \
    --o-pcoa new-TE-soil-coremetrics/locationC_weighted_pcoa.qza 



# location D weighted pcoa
qiime diversity pcoa \
    --i-distance-matrix new-TE-soil-coremetrics/locationD_weighted_distance_matrix.qza \
    --o-pcoa new-TE-soil-coremetrics/locationD_weighted_pcoa.qza 


echo "unweighted visualizations"


# location A emperor unweighted
qiime emperor plot \
    --i-pcoa new-TE-soil-coremetrics/locationA_unweighted_pcoa.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationA_unweighted_unifrac_emperor.qzv


# location B emperor unweighted
qiime emperor plot \
    --i-pcoa new-TE-soil-coremetrics/locationB_unweighted_pcoa.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationB_unweighted_unifrac_emperor.qzv


# location C emperor unweighted
qiime emperor plot \
    --i-pcoa new-TE-soil-coremetrics/locationC_unweighted_pcoa.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationC_unweighted_unifrac_emperor.qzv


# location D emperor unweighted
qiime emperor plot \
    --i-pcoa new-TE-soil-coremetrics/locationD_unweighted_pcoa.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationD_unweighted_unifrac_emperor.qzv


echo "weighted visualizations"

# location A emperor weighted
qiime emperor plot \
    --i-pcoa new-TE-soil-coremetrics/locationA_weighted_pcoa.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationA_weighted_unifrac_emperor.qzv


# location B emperor weighted
qiime emperor plot \
    --i-pcoa new-TE-soil-coremetrics/locationB_weighted_pcoa.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationB_weighted_unifrac_emperor.qzv


# location C emperor weighted
qiime emperor plot \
    --i-pcoa new-TE-soil-coremetrics/locationC_weighted_pcoa.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationC_weighted_unifrac_emperor.qzv


# location D emperor weighted
qiime emperor plot \
    --i-pcoa new-TE-soil-coremetrics/locationD_weighted_pcoa.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationD_weighted_unifrac_emperor.qzv


echo "start unweighted stats"

# location A unweighted type*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/locationA_unweighted_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_type*treatment" \
    --o-visualization new-TE-soil-coremetrics/locationA_unweighted_type_treatment_adonis.qzv 

# location B unweighted type*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/locationB_unweighted_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_type*treatment" \
    --o-visualization new-TE-soil-coremetrics/locationB_unweighted_type_treatment_adonis.qzv 


# location C unweighted type*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/locationC_unweighted_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_type*treatment" \
    --o-visualization new-TE-soil-coremetrics/locationC_unweighted_type_treatment_adonis.qzv 


# location D unweighted type*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/locationD_unweighted_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_type*treatment" \
    --o-visualization new-TE-soil-coremetrics/locationD_unweighted_type_treatment_adonis.qzv 


echo "start weighted stats"

# location A weighted type*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/locationA_weighted_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_type*treatment" \
    --o-visualization new-TE-soil-coremetrics/locationA_weighted_type_treatment_adonis.qzv


# location B weighted type*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/locationB_weighted_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_type*treatment" \
    --o-visualization new-TE-soil-coremetrics/locationB_weighted_type_treatment_adonis.qzv


# location C weighted type*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/locationC_weighted_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_type*treatment" \
    --o-visualization new-TE-soil-coremetrics/locationC_weighted_type_treatment_adonis.qzv


# location D weighted type*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/locationD_weighted_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_type*treatment" \
    --o-visualization new-TE-soil-coremetrics/locationD_weighted_type_treatment_adonis.qzv


echo "filter table by locations"

# location A filter
qiime feature-table filter-samples \
    --i-table new-TE-soil-coremetrics/rarefied_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_A'" \
    --o-filtered-table new-TE-soil-coremetrics/locationA_table.qza


# location B filter
qiime feature-table filter-samples \
    --i-table new-TE-soil-coremetrics/rarefied_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_B'" \
    --o-filtered-table new-TE-soil-coremetrics/locationB_table.qza


# location C filter
qiime feature-table filter-samples \
    --i-table new-TE-soil-coremetrics/rarefied_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_C'" \
    --o-filtered-table new-TE-soil-coremetrics/locationC_table.qza


# location D filter
qiime feature-table filter-samples \
    --i-table new-TE-soil-coremetrics/rarefied_table.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_D'" \
    --o-filtered-table new-TE-soil-coremetrics/locationD_table.qza




##### Filter TE touch by hand

echo "start TE hand beta diversity"

### Filter unweighted DM by hand
qiime diversity filter-distance-matrix \
    --i-distance-matrix new-TE-touch-coremetrics/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_type]= 'skin'" \
    --o-filtered-distance-matrix new-TE-touch-coremetrics/TE_skin_unweighted_distance_matrix.qza



# Filter unweighted DM by PVC
qiime diversity filter-distance-matrix \
    --i-distance-matrix new-TE-touch-coremetrics/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_type] = 'pvc'" \
    --o-filtered-distance-matrix new-TE-touch-coremetrics/TE_pvc_unweighted_distance_matrix.qza



### Filter weighted DM by hand
qiime diversity filter-distance-matrix \
    --i-distance-matrix new-TE-touch-coremetrics/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_type]='skin'" \
    --o-filtered-distance-matrix new-TE-touch-coremetrics/TE_skin_weighted_distance_matrix.qza



# Filter weighted DM by PVC
qiime diversity filter-distance-matrix \
    --i-distance-matrix new-TE-touch-coremetrics/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_type] = 'pvc'" \
    --o-filtered-distance-matrix new-TE-touch-coremetrics/TE_pvc_weighted_distance_matrix.qza



#### Make PCoAs

echo " start unweighted pcoa"

# TE skin unweighted pcoa
qiime diversity pcoa \
    --i-distance-matrix new-TE-touch-coremetrics/TE_skin_unweighted_distance_matrix.qza \
    --o-pcoa new-TE-touch-coremetrics/TE_skin_unweighted_pcoa.qza


# TE pvc unweighted pcoa
qiime diversity pcoa \
    --i-distance-matrix new-TE-touch-coremetrics/TE_pvc_unweighted_distance_matrix.qza \
    --o-pcoa new-TE-touch-coremetrics/TE_pvc_unweighted_pcoa.qza


# TE skin weighted pcoa
qiime diversity pcoa \
    --i-distance-matrix new-TE-touch-coremetrics/TE_skin_weighted_distance_matrix.qza \
    --o-pcoa new-TE-touch-coremetrics/TE_skin_weighted_pcoa.qza


# TE pvc unweighted pcoa
qiime diversity pcoa \
    --i-distance-matrix new-TE-touch-coremetrics/TE_pvc_weighted_distance_matrix.qza \
    --o-pcoa new-TE-touch-coremetrics/TE_pvc_weighted_pcoa.qza



### Make emperor plots

echo "start emperor plots"


# TE skin emperor unweighted
qiime emperor plot \
    --i-pcoa new-TE-touch-coremetrics/TE_skin_unweighted_pcoa.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-touch-coremetrics/TE_skin_unweighted_unifrac_emperor.qzv


# TE pvc emperor unweighted
qiime emperor plot \
    --i-pcoa new-TE-touch-coremetrics/TE_pvc_unweighted_pcoa.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-touch-coremetrics/TE_pvc_unweighted_unifrac_emperor.qzv



# TE skin emperor weighted
qiime emperor plot \
    --i-pcoa new-TE-touch-coremetrics/TE_skin_weighted_pcoa.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-touch-coremetrics/TE_skin_weighted_unifrac_emperor.qzv


# TE pvc emperor weighted
qiime emperor plot \
    --i-pcoa new-TE-touch-coremetrics/TE_pvc_weighted_pcoa.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-touch-coremetrics/TE_pvc_weighted_unifrac_emperor.qzv