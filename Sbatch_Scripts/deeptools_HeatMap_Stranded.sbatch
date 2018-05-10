#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=50G
#SBATCH --job-name="DT_HeatMap_ST"
#SBATCH --output=slurm_%j_DT_HeatMap_ST.out
#SBATCH --mail-user=rluis@medicina.ulisboa.pt
#SBATCH --mail-type=ALL
#SBATCH --image=mcfonsecalab/htstools_plus:latest

cd /home

srun python /home/rluis/MTScripts_DEV/deeptools_HeatMap_Stranded.py $1 --cpu 40 --zMax 2.2