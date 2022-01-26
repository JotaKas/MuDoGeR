#!/bin/bash

# inputs
# 1 is the viral multi fasta file
# 2 is the bins folder
# 3 is the output folder
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database


conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/wish_env


uvigs_file="$1"
bins_folder="$2"
output_folder="$3"

#commands
#1 create output folders
mkdir -p "$output_folder"   
mkdir -p "$output_folder"/potential_host_genomes
mkdir -p "$output_folder"/viral_particles
mkdir -p "$output_folder"/output_results
mkdir -p "$output_folder"/nullmodels

#2 cp and prepare data
cp "$bins_folder"/*fa "$output_folder"/potential_host_genomes
python3 /home/centos/mudoger/MuDoGeR/tools/split-all-seq.py "$uvigs_file" "$output_folder"/viral_particles/viral-particle

#3 build viral model
WIsH -c build -g "$output_folder"/potential_host_genomes/ -m "$output_folder"/modelDir

#4 build bacterial null model
WIsH -c predict -g /mnt/tools/miniconda2/envs/wish-env/database/phages/ -m "$output_folder"/modelDir -r "$output_folder"/nullmodels -b "$output_folder"/nullmodels

#5 run R script to get nullparameters.tsv
cp MuDoGeR/tools/computeNullParameters.R  "$output_folder"/nullmodels
cd "$output_folder"/nullmodels
Rscript computeNullParameters.R 
cd -

#6 run prediction with null model
WIsH -t 20 -c predict -g "$output_folder"/viral_particles/ -m  "$output_folder"/modelDir -r "$output_folder"/output_results/ -b 1 -n "$output_folder"/nullmodels/nullParameters.tsv


