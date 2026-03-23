#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -A ACF-UTK0011
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --ntasks-per-node=12
#SBATCH --job-name=TE-coremetrics
#SBATCH --output=jobout-TE-coremetrics.o%j.txt
#SBATCH --chdir=/lustre/isaac24/proj/UTK0406/EmilyCantrell
#SBATCH --mail-user=sbn572@vols.utk.edu
#SBATCH --mail-type=ALL


# activate conda 
module load anaconda3/2024.06
source activate /lustre/isaac24/proj/UTK0406/shared_conda_envs/qiime2-amplicon-2025.4

THREADS=12

echo "Number of threads to use: $THREADS"

cd /lustre/isaac24/proj/UTK0406/EmilyCantrell



# core metrics - 7k
qiime diversity core-metrics-phylogenetic \
    --i-table taxonomy/SR1_SR2_tax_filter.qza \
    --i-phylogeny tree/tree.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-sampling-depth 7000 \
    --output-dir new-core-metrics-7000


echo "start TE touch"

# TE touch
qiime diversity core-metrics-phylogenetic \
    --i-table merged/TE_touch_table.qza \
    --i-phylogeny tree/tree.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-sampling-depth 10000 \
    --output-dir new-TE-touch-coremetrics



# # Unweighted TE touch adonis stats sample location*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-touch-coremetrics/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*treatment" \
    --o-visualization new-TE-touch-coremetrics/TE_touch_location_treatment_unweighted_adonis.qzv


# # weighted TE touch adonis stats sample location*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-touch-coremetrics/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*treatment" \
    --o-visualization new-TE-touch-coremetrics/TE_touch_location_treatment_weighted_adonis.qzv



# # Unweighted TE touch adonis stats sample location*sample_type
qiime diversity adonis \
    --i-distance-matrix new-TE-touch-coremetrics/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*sample_type" \
    --o-visualization new-TE-touch-coremetrics/TE_touch_location_type_unweighted_adonis.qzv


# # weighted TE touch adonis stats sample location*sample_type
qiime diversity adonis \
    --i-distance-matrix new-TE-touch-coremetrics/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*sample_type" \
    --o-visualization new-TE-touch-coremetrics/TE_touch_location_type_weighted_adonis.qzv


# Unweighted TE pvc adonis stats sample location*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-touch-coremetrics/TE_pvc_unweighted_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*treatment" \
    --o-visualization new-TE-touch-coremetrics/TE_pvc_unweighted_adonis.qzv


# Weighted TE pvc adonis stats sample location*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-touch-coremetrics/TE_pvc_weighted_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*treatment" \
    --o-visualization new-TE-touch-coremetrics/TE_pvc_weighted_adonis.qzv


# Unweighted TE skin adonis stats sample location*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-touch-coremetrics/TE_skin_unweighted_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*treatment" \
    --o-visualization new-TE-touch-coremetrics/TE_skin_unweighted_adonis.qzv


# Weighted TE skin adonis stats sample location*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-touch-coremetrics/TE_skin_weighted_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*treatment" \
    --o-visualization new-TE-touch-coremetrics/TE_skin_weighted_adonis.qzv



echo "Start TE soil"

# TE soil
qiime diversity core-metrics-phylogenetic \
    --i-table merged/TE_soil_table.qza \
    --i-phylogeny tree/tree.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-sampling-depth 27000 \
    --output-dir new-TE-soil-coremetrics


# unweighted Te soil adonis stats - location only
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location" \
    --o-visualization new-TE-soil-coremetrics/TE_soil_location_unweighted_adonis.qzv


# weighted Te soil adonis stats - location only
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location" \
    --o-visualization new-TE-soil-coremetrics/TE_soil_location_weighted_adonis.qzv



# # Unweighted TE soil  adonis stats sample location*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*treatment" \
    --o-visualization new-TE-soil-coremetrics/TE_soil_location_treatment_unweighted_adonis.qzv


# # weighted TE soil adonis stats sample location*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*treatment" \
    --o-visualization new-TE-soil-coremetrics/TE_soil_location_treatment_weighted_adonis.qzv


# Unweighted TE soil  adonis stats sample location*sample_type
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*sample_type" \
    --o-visualization new-TE-soil-coremetrics/TE_soil_location_type_unweighted_adonis.qzv


# weighted TE soil adonis stats sample location*sample_type
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*sample_type" \
    --o-visualization new-TE-soil-coremetrics/TE_soil_location_type_weighted_adonis.qzv

# Unweighted TE soil  adonis stats sample type*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_type*treatment" \
    --o-visualization new-TE-soil-coremetrics/TE_soil_type_treatment_unweighted_adonis.qzv


# weighted TE soil adonis stats sample type*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_type*treatment" \
    --o-visualization new-TE-soil-coremetrics/TE_soil_type_treatment_weighted_adonis.qzv



# unweighted TE soil adonis stats sample type*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*sample_type*treatment" \
    --o-visualization new-TE-soil-coremetrics/TE_soil_location_type_treatment_weighted_adonis.qzv

# weighted TE soil adonis stats sample type*treatment
qiime diversity adonis \
    --i-distance-matrix new-TE-soil-coremetrics/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-formula "sample_location*sample_type*treatment" \
    --o-visualization new-TE-soil-coremetrics/TE_soil_location_type_treatment_weighted_adonis.qzv





######### START ALPHA DIVERSITY

##### EVENNESS

echo "start evenness"

# filter location A 
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/evenness_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_A'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationA_evenness.qza 


# filter location B 
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/evenness_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_B'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationB_evenness.qza 


# filter location C
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/evenness_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_C'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationC_evenness.qza 


# filter location D
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/evenness_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_D'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationD_evenness.qza 


# location A evenness stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationA_evenness.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationA_evenness.qzv


# location B evenness stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationB_evenness.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationB_evenness.qzv


# location C evenness stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationC_evenness.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationC_evenness.qzv


# location D evenness stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationD_evenness.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationD_evenness.qzv



##### FAITH'S DIVERSITY

echo "start faith's"


# filter location A 
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/faith_pd_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_A'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationA_faiths.qza 


# filter location B 
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/faith_pd_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_B'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationB_faiths.qza 


# filter location C
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/faith_pd_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_C'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationC_faiths.qza 


# filter location D
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/faith_pd_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_D'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationD_faiths.qza 


# location A faiths pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationA_faiths.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationA_faiths.qzv


# location B faiths pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationB_faiths.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationB_faiths.qzv


# location C faiths pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationC_faiths.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationC_faiths.qzv


# location D faiths pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationD_faiths.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationD_faiths.qzv


##### RICHNESS

 echo "start richness"


# filter location A 
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/observed_features_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_A'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationA_richness.qza 


# filter location B 
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/observed_features_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_B'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationB_richness.qza 


# filter location C
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/observed_features_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_C'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationC_richness.qza 


# filter location D
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/observed_features_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_D'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationD_richness.qza 


# location A richness pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationA_richness.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationA_richness.qzv


# location B richness pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationB_richness.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationB_richness.qzv


# location C richness pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationC_richness.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationC_richness.qzv


# location D richness pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationD_richness.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationD_richness.qzv



#### SHANNON'S DIVERSITY

 echo "start Shannon's"


# filter location A 
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/shannon_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_A'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationA_shannon.qza 


# filter location B 
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/shannon_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_B'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationB_shannon.qza 


# filter location C
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/shannon_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_C'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationC_shannon.qza 


# filter location D
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-soil-coremetrics/shannon_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_location]='location_D'" \
    --o-filtered-alpha-diversity new-TE-soil-coremetrics/locationD_shannon.qza 


# location A richness pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationA_shannon.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationA_shannon.qzv


# location B richness pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationB_shannon.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationB_shannon.qzv


# location C richness pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationC_shannon.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationC_shannon.qzv


# location D richness pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/locationD_shannon.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/locationD_shannon.qzv



echo "overall AD visualizations"

# overall evenness visualization
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/evenness_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/evenness_vector.qzv


# overall faiths visualization
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/faith_pd_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/faith_pd_vector.qzv


# overall richness visualization
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/observed_features_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/observed_features_vector.qzv



# overall Shannon visualization
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-soil-coremetrics/shannon_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-soil-coremetrics/shannon_vector.qzv




####### TRACE EVIDENCE TOUCH ALPHA DIVERSITY

echo "TE touch Alpha Diversity"

#### Overall TE alpha diversity

# overall evenness visualization
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-touch-coremetrics/evenness_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-touch-coremetrics/evenness_vector.qzv


# overall faiths visualization
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-touch-coremetrics/faith_pd_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-touch-coremetrics/faith_pd_vector.qzv


# overall richness visualization
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-touch-coremetrics/observed_features_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-touch-coremetrics/observed_features_vector.qzv


## TE HAND ALPHA DIVERSITY


# filter hand - faiths
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-touch-coremetrics/faith_pd_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_type]= 'skin'" \
    --o-filtered-alpha-diversity new-TE-touch-coremetrics/TE_hand_faiths.qza 


# filter hand - evenness
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-touch-coremetrics/evenness_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_type]= 'skin'" \
    --o-filtered-alpha-diversity new-TE-touch-coremetrics/TE_hand_evenness.qza 


# filter hand richness
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-touch-coremetrics/observed_features_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_type]= 'skin'" \
    --o-filtered-alpha-diversity new-TE-touch-coremetrics/TE_hand_observedfeatures.qza 


# Hand faiths pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-touch-coremetrics/TE_hand_faiths.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-touch-coremetrics/TE_hand_faiths.qzv

# Hand evenness stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-touch-coremetrics/TE_hand_evenness.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-touch-coremetrics/TE_hand_evenness.qzv


# Hand observed features stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-touch-coremetrics/TE_hand_observedfeatures.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-touch-coremetrics/TE_hand_observedfeatures.qzv






## TE PVC ALPHA DIVERSITY

# filter pvc - faiths
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-touch-coremetrics/faith_pd_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_type]= 'pvc'" \
    --o-filtered-alpha-diversity new-TE-touch-coremetrics/TE_pvc_faiths.qza 


# filter pvc - evenness
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-touch-coremetrics/evenness_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_type]= 'pvc'" \
    --o-filtered-alpha-diversity new-TE-touch-coremetrics/TE_pvc_evenness.qza 


# filter pvc richness
qiime diversity filter-alpha-diversity \
    --i-alpha-diversity new-TE-touch-coremetrics/observed_features_vector.qza \
    --m-metadata-file NIJ_metadata.txt \
    --p-where "[sample_type]= 'pvc'" \
    --o-filtered-alpha-diversity new-TE-touch-coremetrics/TE_pvc_observedfeatures.qza 


# pvc faiths pd stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-touch-coremetrics/TE_pvc_faiths.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-touch-coremetrics/TE_pvc_faiths.qzv

# pvc evenness stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-touch-coremetrics/TE_pvc_evenness.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-touch-coremetrics/TE_pvc_evenness.qzv


# pvc observed features stats
qiime diversity alpha-group-significance \
    --i-alpha-diversity new-TE-touch-coremetrics/TE_pvc_observedfeatures.qza \
    --m-metadata-file NIJ_metadata.txt \
    --o-visualization new-TE-touch-coremetrics/TE_pvc_observedfeatures.qzv
