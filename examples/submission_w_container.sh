#!/bin/bash
#SBATCH --job-name=dev
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=30GB 

module load apptainer/1.2.4

# path to snakefile and config 
snakefile=</path/to/snakefile>
config=</path/to/configfile>

# path to container and working directory
image=</path/to/container>
working_dir=</path/to/working_directory>

# unlock directory
apptainer exec --contain --cleanenv --pwd "$working_dir" $image snakemake -s ${snakefile} --configfile ${config} --unlock 

# Run pipeline
apptainer exec --contain --cleanenv --pwd "$working_dir" $image snakemake -s ${snakefile} --configfile ${config} --cores 5