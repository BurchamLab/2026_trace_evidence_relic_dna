#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -A ACF-UTK0011
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --ntasks-per-node=12
#SBATCH --job-name=DADA2
#SBATCH --output=jobout-DADA2.o%j.txt
#SBATCH --chdir=/lustre/isaac24/proj/UTK0406/EmilyCantrell
#SBATCH --mail-user=sbn572@vols.utk.edu
#SBATCH --mail-type=ALL


# activate conda 
module load anaconda3/2024.06
source activate /lustre/isaac24/proj/UTK0406/shared_conda_envs/qiime2-amplicon-2025.7

THREADS=12

echo "Number of threads to use: $THREADS"

cd /lustre/isaac24/proj/UTK0406/EmilyCantrell

mkdir -p /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR1_denoise
mkdir -p /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR2_denoise


 # DADA2 paired SR1 denoising 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs import/SR1_demux.qza \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 240 \
  --p-n-threads $THREADS \
  --o-table /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR1_denoise/dada2_table_pe240.qza \
  --o-representative-sequences /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR1_denoise/SR1_dada2_rep_set_pe240.qza \
  --o-denoising-stats /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR1_denoise/SR1_dada2_stats_pe240.qza
  
qiime metadata tabulate \
  --m-input-file /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR1_denoise/SR1_dada2_stats_pe240.qza  \
  --o-visualization /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR1_denoise/SR1_dada2_stats_pe240.qzv


qiime feature-table summarize \
  --i-table /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR1_denoise/dada2_table_pe240.qza \
  --o-visualization /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR1_denoise/SR1_dada2_table_pe240.qzv


#DADA2 paired SR2 denoising
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs import/SR2_demux.qza \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 240 \
  --p-n-threads $THREADS \
  --o-table /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR2_denoise/SR2_dada2_table_pe240.qza \
  --o-representative-sequences /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR2_denoise/SR2_dada2_rep_set_pe240.qza \
  --o-denoising-stats /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR2_denoise/SR2_dada2_stats_pe240.qza
  
qiime metadata tabulate \
  --m-input-file /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR2_denoise/SR2_dada2_stats_pe240.qza  \
  --o-visualization /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR2_denoise/SR2_dada2_stats_pe240.qzv
    

qiime feature-table summarize \
  --i-table /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR2_denoise/SR2_dada2_table_pe240.qza \
  --o-visualization /lustre/isaac24/proj/UTK0406/EmilyCantrell/SR2_denoise/SR2_dada2_table_pe240.qzv



# SR1 single end denoise
qiime dada2 denoise-single \
  --i-demultiplexed-seqs /lustre/isaac24/proj/UTK0406/EmilyCantrell/import/SR1_single.qza \
  --p-trunc-len 240 \
  --p-n-threads $THREADS \
  --o-representative-sequences SR1_denoise/SR1_dada2_single_repseqs.qza \
  --o-table SR1_denoise/SR1_dada2_single_table.qza \
  --o-denoising-stats SR1_denoise/SR1_dada2_single_stats.qza

# stats visualization
qiime metadata tabulate \
  --m-input-file SR1_denoise/SR1_dada2_single_stats.qza \
  --o-visualization SR1_denoise/SR1_dada2_single_stats.qzv

# table visualization
qiime feature-table summarize \
  --i-table SR1_denoise/SR1_dada2_single_table.qza \
  --o-visualization SR1_denoise/SR1_dada2_single_table.qzv


# rep seqs visualization
qiime feature-table tabulate-seqs \
  --i-data SR1_denoise/SR1_dada2_single_repseqs.qza \
  --o-visualization SR1_denoise/SR1_dada2_single_repseqs.qzv 


#SR2 single end denoise
qiime dada2 denoise-single \
  --i-demultiplexed-seqs /lustre/isaac24/proj/UTK0406/EmilyCantrell/import/SR2_single.qza \
  --p-trunc-len 240 \
  --p-n-threads $THREADS \
  --o-representative-sequences SR2_denoise/SR2_dada2_single_repseqs.qza \
  --o-table SR2_denoise/SR2_dada2_single_table.qza \
  --o-denoising-stats SR2_denoise/SR2_dada2_single_stats.qza


# stats visualization
qiime metadata tabulate \
  --m-input-file SR2_denoise/SR2_dada2_single_stats.qza \
  --o-visualization SR2_denoise/SR2_dada2_single_stats.qzv


# table visualization
qiime feature-table summarize \
  --i-table SR2_denoise/SR2_dada2_single_table.qza \
  --o-visualization SR2_denoise/SR2_dada2_single_table.qzv


# rep seqs visualization
qiime feature-table tabulate-seqs \
  --i-data SR2_denoise/SR2_dada2_single_repseqs.qza \
  --o-visualization SR2_denoise/SR2_dada2_single_repseqs.qzv 