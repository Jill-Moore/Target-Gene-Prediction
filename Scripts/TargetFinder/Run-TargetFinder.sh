#!/bin/bash
#SBATCH --nodes 1
#SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH --output=/home/moorej3/Job-Logs/jobid_%A.output
#SBATCH --error=/home/moorej3/Job-Logs/jobid_%A.error
#SBATCH --partition=12hours

data=NewCalls-HiC
version=v1

setDir=~/Lab/Target-Gene/Benchmark
train=$setDir/$data-Benchmark.$version.txt
peakDir=/data/zusers/garnickk/targetfinder/GM12878/peaks/
featureDir=~/Lab/Target-Gene/Target-Finder/Feature-Matrices
outputDir=~/Lab/Target-Gene/Target-Finder/Results
enhancers=~/Lab/Target-Gene/Target-Finder/GM12878-Enhancers.bed
scriptDir=~/Projects/Target-Gene-Prediction/Scripts/TargetFinder
tss=~/Lab/Reference/Human/hg19/Gencode19/TSS.Filtered.bed
bedtools=~/bin/bedtools2/bin/bedtools

mkdir -p /tmp/moorej3/$SLURM_JOBID-$jid
cd /tmp/moorej3/$SLURM_JOBID-$jid

mkdir -p $outputDir

######## Creating Enhancer Feature Matrix ################
if [ ! -f "$featureDir/$data-Enhancer-Feature-Matrix.txt" ]
then
    echo -e "Generating enhancer feature matrix..."
    cat $train | awk '{print $1}' | sort -u  > cres
    awk 'FNR==NR {x[$1];next} ($4 in x)' cres $enhancers | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $9 }' > enhancers
    for k in $(seq 99)
    do
        echo $k
        peakFile=$(cat ~/Lab/Target-Gene/Target-Finder/Dataset-Peak-List.txt | awk -F "\t" '{if (NR == '$k') \
            print $1}')
    
        $bedtools intersect -wo -a enhancers -b $peakDir/$peakFile > tmp
        python $scriptDir/process.overlaps.py enhancers tmp | sort -k1,1 | \
            awk 'BEGIN {print "cREs" "\t" "'$peakFile'"}{print $0}'> col.$k
    done
    paste col.* | awk '{printf "%s\t", $1; for(i=2;i<=NF;i+=2) printf "%s\t", \
        $i;print ""}' > $featureDir/$data-Enhancer-Feature-Matrix.txt
    rm col.*
fi

######## Creating TSS Feature Matrix ################
if [ ! -f "$featureDir/$data-TSS-Feature-Matrix.txt" ]
then
    echo -e "Generating tss feature matrix..."
    cat $train | awk '{print $2}' | sort -u  > genes
    awk 'FNR==NR {x[$1];next} ($7 in x)' genes $tss | \
        awk '{print $1 "\t" $2-500 "\t" $3+500 "\t" $4 "\t" $7 }' > tss
    for k in $(seq 99)
    do
        echo $k
        peakFile=$(cat ~/Lab/Target-Gene/Target-Finder/Dataset-Peak-List.txt | awk -F "\t" '{if (NR == '$k') \
             print $1}')
        $bedtools intersect -wo -a tss -b $peakDir/$peakFile > tmp
        python $scriptDir/process.overlaps.py tss tmp | sort -k1,1 | \
            awk 'BEGIN {print "cREs" "\t" "'$peakFile'"}{print $0}'> col.$k
    done
    paste col.* | awk '{printf "%s\t", $1; for(i=2;i<=NF;i+=2) printf "%s\t",\
        $i;print ""}' > $featureDir/$data-TSS-Feature-Matrix.txt
    rm col.*
fi

######## Creating Window Matrix ################
if [ ! -f "$featureDir/$data-Window-Feature-Matrix.txt" ]
then
    echo -e "Generating window feature matrix..."
    cat $train | sort -u > pairs
    python $scriptDir/create.window.py $tss $enhancers pairs > windows
    for k in $(seq 99)
    do
        echo $k
        peakFile=$(cat ~/Lab/Target-Gene/Target-Finder/Dataset-Peak-List.txt | awk -F "\t" '{if (NR == '$k') \
             print $1}')
        $bedtools intersect -wo -a windows -b $peakDir/$peakFile > tmp
        python $scriptDir/process.overlaps.py windows tmp | sort -k1,1 | \
            awk 'BEGIN {print "cREs" "\t" "'$peakFile'"}{print $0}'> col.$k
    done
    paste col.* | awk '{printf "%s\t", $1; for(i=2;i<=NF;i+=2) printf "%s\t",\
        $i;print ""}' > $featureDir/$data-Window-Feature-Matrix.txt
    rm col.*
fi

######## Creating Distance Matrix ################
if [ ! -f "$featureDir/$data-Distance.txt" ]
then
    echo -e "Generating distance matrix..."
    cat $train | sort -u > pairs
    python $scriptDir/calculate.distance.py $tss $enhancers pairs > \
        $featureDir/$data-Distance.txt 
fi


######## Running Random Forest/GBM ################
echo -e "Running Model..."
cat $train | awk '{print $2}' | sort -u  > genes
awk 'FNR==NR {x[$1];next} ($7 in x)' genes $tss | \
awk '{print $1 "\t" $2-500 "\t" $3+500 "\t" $4 "\t" $7 }' > tss

python $scriptDir/random.forest.py $train $featureDir/$data-Enhancer-Feature-Matrix.txt \
    $featureDir/$data-TSS-Feature-Matrix.txt $featureDir/$data-Window-Feature-Matrix.txt \
    tss $data $outputDir $version

#python gbm.targetfinder.window.py $train $val $featureDir/$data-Enhancer-Feature-Matrix.txt \
#    $featureDir/$data-TSS-Feature-Matrix.txt $featureDir/$data-Window-Feature-Matrix.txt \
#    tss $data > $data-GBM-Features.txt

rm -r /tmp/moorej3/$SLURM_JOBID-$jid
