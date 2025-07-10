#!/bin/bash

############# SLURM SETTINGS #############
#SBATCH --account=none
#SBATCH --job-name=run_dattr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=nodes
#SBATCH --time=0-10:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=1

############# LOADING MODULES (optional) #############
module load apps/R
module load apps/R/4.4.1/gcc-8.5.0+openblas-0.3.28
export LD_LIBRARY_PATH=/users/2850244m/Documents/install_multinma/condaenv/lib:${LD_LIBRARY_PATH}

############# MY CODE #############
Rscript obj2.R
