#!/bin/bash
#SBATCH --job-name=scXpress
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=30GB 

# path to snakefile and config 
snakefile=</path/to/snakefile>
config=</path/to/configfile>

# unlock directory 
 snakemake -s ${snakefile} --configfile ${config} --unlock 

# run pipeline
 snakemake -s ${snakefile} --configfile ${config} --cores 5 