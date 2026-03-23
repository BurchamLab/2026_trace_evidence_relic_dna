#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -A ACF-UTK0011
#SBATCH --partition=campus-bigmem
#SBATCH --qos=campus-bigmem
#SBATCH --ntasks-per-node=12
#SBATCH --job-name=Taxonomy
#SBATCH --output=jobout-taxonomy.o%j.txt
#SBATCH --chdir=/lustre/isaac24/proj/UTK0406/EmilyCantrell
#SBATCH --mail-user=sbn572@vols.utk.edu
#SBATCH --mail-type=ALL


# activate conda 
module load anaconda3/2024.06
source activate /lustre/isaac24/proj/UTK0406/shared_conda_envs/qiime2-amplicon-2025.7

THREADS=12

echo "Number of threads to use: $THREADS"


cd /lustre/isaac24/proj/UTK0406/EmilyCantrell

# download SEPP reference
mkdir -p /lustre/isaac24/proj/UTK0406/EmilyCantrell/tree

wget -O "tree/sepp-refs-silva-128.qza" "https://data.qiime2.org/classifiers/sepp-ref-dbs/sepp-refs-silva-128.qza"

# generate tree
qiime fragment-insertion sepp \
  --i-representative-sequences merged/SR1_SR2_merged_rep_seqs.qza \
  --i-reference-database tree/sepp-refs-silva-128.qza \
  --o-tree tree/tree.qza \
  --o-placements tree/tree_placements.qza \
  --p-threads $THREADS
  
# download greengenes2 515f/806r classifier
mkdir -p /lustre/isaac24/proj/UTK0406/EmilyCantrell/taxonomy

wget -O "taxonomy/2024.09.backbone.v4.nb.sklearn-1.4.2.qza" "https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes2/2024.09.backbone.v4.nb.sklearn-1.4.2.qza"

# classify sequences
qiime feature-classifier classify-sklearn \
  --i-reads merged/SR1_SR2_merged_rep_seqs.qza \
  --i-classifier taxonomy/2024.09.backbone.v4.nb.sklearn-1.4.2.qza \
  --p-n-jobs $THREADS \
  --o-classification taxonomy/taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy/taxonomy.qza \
  --o-visualization taxonomy/taxonomy.qzv
  

qiime taxa barplot \
  --i-table merged/SR1_SR2_merged_table.qza \
  --i-taxonomy taxonomy/taxonomy.qza \
  --m-metadata-file NIJ_metadata.txt \
  --o-visualization taxonomy/taxa-bar-plots.qzv


# filter based on taxonomy
qiime taxa filter-table \
  --i-table merged/SR1_SR2_merged_table.qza \
  --i-taxonomy taxonomy/taxonomy.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table taxonomy/SR1_SR2_tax_filter.qza

qiime feature-table summarize \
	--i-table taxonomy/SR1_SR2_tax_filter.qza \
 	--m-sample-metadata-file NIJ_metadata.txt \
 	--o-visualization taxonomy/SR1_SR2_tax_filter.qzv