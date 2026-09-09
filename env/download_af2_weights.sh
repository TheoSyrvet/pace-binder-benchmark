#!/bin/bash
#SBATCH --job-name=af2_weights
#SBATCH --partition=shared-cpu
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --output=%x_%j.log

mkdir -p ~/master_thesis/BindCraft/repo/params
cd ~/master_thesis/BindCraft/repo/params
wget https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
tar -xvf alphafold_params_2022-12-06.tar
rm alphafold_params_2022-12-06.tar
