#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH -A ACF-UTK0011
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --ntasks-per-node=12
#SBATCH --job-name=import
#SBATCH --output=jobout-import.o%j.txt
#SBATCH --chdir=/lustre/isaac24/proj/UTK0406/EmilyCantrell
#SBATCH --mail-user=sbn572@vols.utk.edu
#SBATCH --mail-type=ALL


# activate conda 
module load anaconda3/2024.06
source activate /lustre/isaac24/proj/UTK0406/shared_conda_envs/qiime2-amplicon-2025.7

THREADS=12

echo "Number of threads to use: $THREADS"

cd /lustre/isaac24/proj/UTK0406/EmilyCantrell

mkdir -p /lustre/isaac24/proj/UTK0406/EmilyCantrell/import

# export SR1
qiime tools export \
    --input-path /lustre/isaac24/proj/UTK0406/raw_data/nij/import/SR1_merged.qza \
    --output-path import/SR1_exported


# export SR2
qiime tools export \
    --input-path /lustre/isaac24/proj/UTK0406/raw_data/nij/import/SR2_merged.qza \
    --output-path import/SR2_exported


# SR1 import
qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path import/SR1_manifest.txt \
    --output-path import/SR1_single.qza \
    --input-format SingleEndFastqManifestPhred33V2



# SR2 import
qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path import/SR2_manifest.txt \
    --output-path import/SR2_single.qza \
    --input-format SingleEndFastqManifestPhred33V2