#!/bin/bash
#SBATCH --nodes 1
#SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH --output=/home/moorej3/Job-Logs/jobid_%A.output
#SBATCH --error=/home/moorej3/Job-Logs/jobid_%A.error
#SBATCH --partition=12hours

data=Ruan-CTCF
version=v1

setDir=~/Lab/Target-Gene/Benchmark
train=$setDir/$data-Benchmark.$version.txt
featureDir=~/Lab/Target-Gene/PEP/Feature-Matrices
outputDir=~/Lab/Target-Gene/Target-Finder/Results
enhancers=~/Lab/Target-Gene/Target-Finder/GM12878-Enhancers.bed
scriptDir=~/Projects/Target-Gene-Prediction/Scripts/PEP
tss=~/Lab/Reference/Human/hg19/Gencode19/TSS.Filtered.bed

genome=~/Lab/Reference/Human/hg19/hg19.2bit
fimo=~/bin/meme/bin/fimo
motifs=~/Lab/Reference/Motifs/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
outputDir=~/Lab/Target-Gene/PEP/Results

mkdir -p /tmp/moorej3/$SLURM_JOBID-$jid
cd /tmp/moorej3/$SLURM_JOBID-$jid

mkdir -p $outputDir

######## Creating Enhancer Feature Matrix ################
if [ ! -f "$featureDir/$data-Enhancer-Feature-Matrix.txt" ]
then
    echo -e "Generating enhancer feature matrix..."
    cat $train | awk '{print $1}' | sort -u  > cres
    awk 'FNR==NR {x[$1];next} ($4 in x)' cres $enhancers | \
    awk '{print $1 "\t" $2-4000 "\t" $3+4000 "\t" $4}' > enhancers
    ~/bin/twoBitToFa $genome enhancers.fa -bed=enhancers
    $fimo --text $motifs enhancers.fa  > fimo.out
    python $scriptDir/count.motifs.py HOCOMOCO-Motif-Header.txt enhancers fimo.out \
    > $featureDir/$data-Enhancer-Feature-Matrix.txt
fi

######## Creating TSS Feature Matrix ################
if [ ! -f "$featureDir/$data-TSS-Feature-Matrix.txt" ]
then
    echo -e "Generating tss feature matrix..."
    cat $train | awk '{print $2}' | sort -u  > genes
    awk 'FNR==NR {x[$1];next} ($7 in x)' genes $tss | \
        awk '{print $1 "\t" $2-500 "\t" $3+500 "\t" $4}' > tss
    ~/bin/twoBitToFa $genome tss.fa -bed=tss
    $fimo --text $motifs tss.fa  > fimo.out
    python $scriptDir/count.motifs.py HOCOMOCO-Motif-Header.txt tss fimo.out \
    > $featureDir/$data-TSS-Feature-Matrix.txt
fi

######## Creating Distance Matrix ################
#if [ ! -f "$featureDir/$data-Distance.txt" ]
#then
#    echo -e "Generating distance matrix..."
#    cat $val $train | sort -u > pairs
#    python calculate.distance.py $tss $enhancers pairs > \
#        $featureDir/$data-Distance.txt 
#fi

######## Running Random Forest ################
echo -e "Running Model..."
cat $train | awk '{print $2}' | sort -u  > genes
awk 'FNR==NR {x[$1];next} ($7 in x)' genes $tss | \
awk '{print $1 "\t" $2-500 "\t" $3+500 "\t" $4 "\t" $7 }' > tss

python $scriptDir/xgb.py $train $featureDir/$data-Enhancer-Feature-Matrix.txt \
    $featureDir/$data-TSS-Feature-Matrix.txt tss $data $outputDir $version 

rm -r /tmp/moorej3/$SLURM_JOBID-$jid
