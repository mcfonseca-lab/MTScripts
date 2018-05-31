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

full_name1=$(readlink -f $1)
cd /home/rluis/scratch

mkdir DT_HeatMap_ST_$SLURM_JOB_ID
cd DT_HeatMap_ST_$SLURM_JOB_ID


cp $full_name1 .

file1=$(basename $1)

srun shifter python /home/rluis/MTScripts_DEV/deeptools_TSS1_TSS2_metaProfile_Stranded.py $file1 --cpu 40 --yMax 20
rm $file1
mv * $SLURM_SUBMIT_DIR
cd ..
rmdir DT_HeatMap_ST_$SLURM_JOB_ID
