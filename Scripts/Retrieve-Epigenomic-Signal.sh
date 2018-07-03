#!/bin/bash
#SBATCH --nodes 1
#SBATCH --time=00:30:00
#SBATCH --mem=5G
#SBATCH --array=1-149%20
#SBATCH --output=/home/moorej3/Job-Logs/jobid_%A_%a.output
#SBATCH --error=/home/moorej3/Job-Logs/jobid_%A_%a.error

j=$SLURM_ARRAY_TASK_ID

dir=~/Lab/Target-Gene/Benchmark/Properties/TF-Enrichment-ELS
dataDir=/data/projects/encode/data/

files=$dir/GM12878-TF-List.05-14-18.txt
ccREs=$dir/All-ELS.bed

exp=$(awk -F "\t" '{if (NR=='$j') print $1}' $files)
sig=$(awk -F "\t" '{if (NR=='$j') print $2}' $files)
target=$(awk -F "\t" '{if (NR=='$j') print $3}' $files)

mkdir -p /tmp/moorej3/$SLURM_JOBID"-"$SLURM_ARRAY_TASK_ID
cd /tmp/moorej3/$SLURM_JOBID"-"$SLURM_ARRAY_TASK_ID

/bin/sleep   `/usr/bin/expr $RANDOM % 60`
~/bin/bigWigAverageOverBed $dataDir/$exp/$sig.bigWig $ccREs out

mv out $dir/signal-output/$exp"-"$sig.tab
rm -r /tmp/moorej3/$SLURM_JOBID"-"$SLURM_ARRAY_TASK_ID

