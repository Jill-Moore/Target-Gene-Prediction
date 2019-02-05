#!/bin/bash
#SBATCH --nodes 1
#SBATCH --time=04:00:00
#SBATCH --mem=10G
#SBATCH --output=/home/moorej3/Job-Logs/jobid_%A.output
#SBATCH --error=/home/moorej3/Job-Logs/jobid_%A.error
#SBATCH --partition=4hours

data1=Ruan-CTCF
#data2=K562.HiC

#datasets=("Aiden-HiC" "HeLa.HiC" "HMEC.HiC" "IMR90.HiC" "K562.HiC" "NHEK.HiC")
datasets=("Ruan-CTCF" "HeLa.Ruan-CTCF")

for data2 in ${datasets[@]}
do
echo $data2

version=v1

setDir=~/Lab/Target-Gene/Benchmark
train=$setDir/$data2-Benchmark.$version.txt
featureDir=~/Lab/Target-Gene/Target-Finder/Feature-Matrices
outputDir1=~/Lab/Target-Gene/Target-Finder/Results-4C
outputDir2=~/Lab/Target-Gene/Target-Finder/Results-CrossCell
scriptDir=~/Projects/Target-Gene-Prediction/Scripts/TargetFinder
tss=~/Lab/Reference/Human/hg19/Gencode19/TSS.Filtered.bed
ccREs=~/Lab/ENCODE/Encyclopedia/V4/Registry/V4-hg19/hg19-ccREs-Simple.bed

mkdir -p /tmp/moorej3/$SLURM_JOBID-$jid
cd /tmp/moorej3/$SLURM_JOBID-$jid

mkdir -p $outputDir2

######## Adjust chrom CV ################
echo -e "Adjusting chrom CV"

train1=$setDir/$data1-Benchmark.$version.txt
cat $train1 | awk '{print $1 "\t" $4}' | sort -u | sort -k1,1  > cres
awk 'FNR==NR {x[$1];next} ($5 in x)' cres $ccREs | \
    awk '{print $1 "\t" $5}' | sort -k2,2 > tmp1

paste tmp1 cres | awk '{print $1 "\t" $4}' | sort -u  > chrom-cv.txt
wc -l chrom-cv.txt
awk 'FNR==NR {x[$1];next} ($5 in x)' $train $ccREs > enhancers

python $scriptDir/reassign.chromCV.py $train enhancers chrom-cv.txt > new-train

######## Running Random Forest/GBM ################
echo -e "Running Model..."
cat new-train | awk '{print $2}' | sort -u  > genes
awk 'FNR==NR {x[$1];next} ($7 in x)' genes $tss | \
awk '{print $1 "\t" $2-500 "\t" $3+500 "\t" $4 "\t" $7 }' > tss

python $scriptDir/gradient.boosting.CC.py new-train $featureDir/$data2-4C-Enhancer-Feature-Matrix.txt \
    $featureDir/$data2-4C-TSS-Feature-Matrix.txt $featureDir/$data2-4C-Window-Feature-Matrix.txt \
    $featureDir/$data2-Distance.txt tss $data1 $data2 $outputDir1 $outputDir2 $version


rm -r /tmp/moorej3/$SLURM_JOBID-$jid/*

done

rm -r /tmp/moorej3/$SLURM_JOBID-$jid/
