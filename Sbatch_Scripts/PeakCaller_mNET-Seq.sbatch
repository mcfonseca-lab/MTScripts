#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=240G
#SBATCH --job-name="PeakCall_NET-Seq"
#SBATCH --output=slurm_%j_PeakCall_NET-Seq.out
#SBATCH --mail-user=rluis@medicina.ulisboa.pt
#SBATCH --mail-type=ALL
#SBATCH --image=mcfonsecalab/htstools_plus:latest

################################
##			      ##
##  Arg1 -> bedFile           ##
##  Arg2 -> bamFile           ##
##  Arg3 -> output File Name  ##
##			      ##
################################


bedFile=$(readlink -f $1)
bamFile=$(readlink -f $2)
cd /home/rluis/scratch

mkdir PeakCall_NET-Seq_$SLURM_JOB_ID
cd PeakCall_NET-Seq_$SLURM_JOB_ID


cp $bedFile .
cp $bamFile .

basenameBedFile=$(basename $bedFile)
basenameBamFile=$(basename $bamFile)

srun shifter python /home/rluis/MTScripts_DEV/PeakCaller_mNET-Seq.py $basenameBedFile $basenameBamFile --outNameFile $3

rm $bedFile
rm $bamFile

mv * $SLURM_SUBMIT_DIR
cd ..
rmdir PeakCall_NET-Seq_$SLURM_JOB_ID

