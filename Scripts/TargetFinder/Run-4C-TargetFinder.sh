#!/bin/bash
#SBATCH --nodes 1
#SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH --output=/home/moorej3/Job-Logs/jobid_%A.output
#SBATCH --error=/home/moorej3/Job-Logs/jobid_%A.error
#SBATCH --partition=12hours

data=Ruan-RNAPII
biosample=GM12878
version=v1

setDir=~/Lab/Target-Gene/Benchmark
train=$setDir/$data-Benchmark.$version.txt
featureDir=~/Lab/Target-Gene/Target-Finder/Feature-Matrices
outputDir=~/Lab/Target-Gene/Target-Finder/Results-4C
ccREs=~/Lab/ENCODE/Encyclopedia/V4/Registry/V4-hg19/hg19-ccREs-Simple.bed
#enhancers=~/Lab/Target-Gene/Target-Finder/GM12878-Enhancers.bed
scriptDir=~/Projects/Target-Gene-Prediction/Scripts/TargetFinder
tss=~/Lab/Reference/Human/hg19/Gencode19/TSS.Filtered.bed
bedtools=~/bin/bedtools2/bin/bedtools
tfList=~/Lab/Target-Gene/Target-Finder/Four-Core-TFs.txt

mkdir -p /tmp/moorej3/$SLURM_JOBID-$jid
cd /tmp/moorej3/$SLURM_JOBID-$jid

mkdir -p $outputDir

######## Creating Enhancer Feature Matrix ################
if [ ! -f "$featureDir/$data-4C-Enhancer-Feature-Matrix.txt" ]
then
    echo -e "Generating enhancer feature matrix..."
    cat $train | awk '{print $1}' | sort -u  > cres
    awk 'FNR==NR {x[$1];next} ($5 in x)' cres $ccREs | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 }' > enhancers
    
    d=$(awk '{if ($9 == "'$biosample'") print $1"/"$2}' $tfList)    
    cp /data/projects/encode/data/$d.bed.gz bed.gz
    gunzip bed.gz
    $bedtools intersect -wo -a enhancers -b bed > tmp
    python $scriptDir/process.overlaps.py enhancers tmp peaks | sort -k1,1 | \
        awk 'BEGIN {print "ccREs" "\t" "DNase"}{print $0}'> col.1
    rm bed    

    hme3=$(awk '{if ($9 == "'$biosample'") print $3"/"$4}' $tfList)
    cp /data/projects/encode/data/$hme3.bed.gz bed.gz
    gunzip bed.gz
    $bedtools intersect -wo -a enhancers -b bed > tmp
    python $scriptDir/process.overlaps.py enhancers tmp peaks | sort -k1,1 | \
        awk 'BEGIN {print "ccREs" "\t" "H3K4me3"}{print $0}'> col.2
    rm bed

    h27ac=$(awk '{if ($9 == "'$biosample'") print $5"/"$6}' $tfList)
    cp /data/projects/encode/data/$h27ac.bed.gz bed.gz
    gunzip bed.gz
    $bedtools intersect -wo -a enhancers -b bed > tmp
    python $scriptDir/process.overlaps.py enhancers tmp peaks | sort -k1,1 | \
        awk 'BEGIN {print "ccREs" "\t" "H3K27ac"}{print $0}'> col.3
    rm bed

    ctcf=$(awk '{if ($9 == "'$biosample'") print $7"/"$8}' $tfList)
    cp /data/projects/encode/data/$ctcf.bed.gz bed.gz
    gunzip bed.gz
    $bedtools intersect -wo -a enhancers -b bed > tmp
    python $scriptDir/process.overlaps.py enhancers tmp peaks | sort -k1,1 | \
        awk 'BEGIN {print "ccREs" "\t" "CTCF"}{print $0}'> col.4
    rm bed 

    paste col.* | awk '{printf "%s\t", $1; for(i=2;i<=NF;i+=2) printf "%s\t", \
        $i;print ""}' > $featureDir/$data-4C-Enhancer-Feature-Matrix.txt
    rm col.*
fi

######## Creating TSS Feature Matrix ################
if [ ! -f "$featureDir/$data-4C-TSS-Feature-Matrix.txt" ]
then
    echo -e "Generating tss feature matrix..."
    cat $train | awk '{print $2}' | sort -u  > genes
    awk 'FNR==NR {x[$1];next} ($7 in x)' genes $tss | \
        awk '{print $1 "\t" $2-500 "\t" $3+500 "\t" $4 "\t" $7 }' > tss
    
    d=$(awk '{if ($9 == "'$biosample'") print $1"/"$2}' $tfList)
    cp /data/projects/encode/data/$d.bed.gz bed.gz
    gunzip bed.gz
    $bedtools intersect -wo -a tss -b bed > tmp
    python $scriptDir/process.overlaps.py tss tmp peaks | sort -k1,1 | \
        awk 'BEGIN {print "cREs" "\t" "DNase"}{print $0}'> col.1
    rm bed

    hme3=$(awk '{if ($9 == "'$biosample'") print $3"/"$4}' $tfList)
    cp /data/projects/encode/data/$hme3.bed.gz bed.gz
    gunzip bed.gz
    $bedtools intersect -wo -a tss -b bed > tmp
    python $scriptDir/process.overlaps.py tss tmp peaks | sort -k1,1 | \
        awk 'BEGIN {print "cREs" "\t" "H3K4me3"}{print $0}'> col.2
    rm bed

    h27ac=$(awk '{if ($9 == "'$biosample'") print $5"/"$6}' $tfList)
    cp /data/projects/encode/data/$h27ac.bed.gz bed.gz
    gunzip bed.gz
    $bedtools intersect -wo -a tss -b bed > tmp
    python $scriptDir/process.overlaps.py tss tmp peaks | sort -k1,1 | \
        awk 'BEGIN {print "cREs" "\t" "H3K27ac"}{print $0}'> col.3
    rm bed

    ctcf=$(awk '{if ($9 == "'$biosample'") print $7"/"$8}' $tfList)
    cp /data/projects/encode/data/$ctcf.bed.gz bed.gz
    gunzip bed.gz
    $bedtools intersect -wo -a tss -b bed > tmp
    python $scriptDir/process.overlaps.py tss tmp peaks | sort -k1,1 | \
        awk 'BEGIN {print "ccREs" "\t" "CTCF"}{print $0}'> col.4
    rm bed
    
    paste col.* | awk '{printf "%s\t", $1; for(i=2;i<=NF;i+=2) printf "%s\t",\
        $i;print ""}' > $featureDir/$data-4C-TSS-Feature-Matrix.txt
    rm col.*
fi

######## Creating Window Matrix ################
if [ ! -f "$featureDir/$data-4C-Window-Feature-Matrix.txt" ]
then
    echo -e "Generating window feature matrix..."
    cat $train | sort -u > pairs
    awk 'FNR==NR {x[$1];next} ($5 in x)' pairs $ccREs | \
        awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 }' > enhancers

    python $scriptDir/create.window.py $tss enhancers pairs > windows
    
    d=$(awk '{if ($9 == "'$biosample'") print $1"/"$2}' $tfList)
    cp /data/projects/encode/data/$d.bed.gz bed.gz
    gunzip bed.gz
    $bedtools intersect -wo -a windows -b bed > tmp
    python $scriptDir/process.overlaps.py windows tmp peaks | sort -k1,1 | \
        awk 'BEGIN {print "cREs" "\t" "DNase"}{print $0}'> col.1
    rm bed

    hme3=$(awk '{if ($9 == "'$biosample'") print $3"/"$4}' $tfList)
    cp /data/projects/encode/data/$hme3.bed.gz bed.gz
    gunzip bed.gz
    $bedtools intersect -wo -a windows -b bed > tmp
    python $scriptDir/process.overlaps.py windows tmp peaks | sort -k1,1 | \
        awk 'BEGIN {print "cREs" "\t" "H3K4me3"}{print $0}'> col.2
    rm bed

    h27ac=$(awk '{if ($9 == "'$biosample'") print $5"/"$6}' $tfList)
    cp /data/projects/encode/data/$h27ac.bed.gz bed.gz
    gunzip bed.gz
    $bedtools intersect -wo -a windows -b bed > tmp
    python $scriptDir/process.overlaps.py windows tmp peaks | sort -k1,1 | \
        awk 'BEGIN {print "cREs" "\t" "H3K27ac"}{print $0}'> col.3
    rm bed

    ctcf=$(awk '{if ($9 == "'$biosample'") print $7"/"$8}' $tfList)
    cp /data/projects/encode/data/$ctcf.bed.gz bed.gz
    gunzip bed.gz
    $bedtools intersect -wo -a windows -b bed > tmp
    python $scriptDir/process.overlaps.py windows tmp peaks | sort -k1,1 | \
        awk 'BEGIN {print "ccREs" "\t" "CTCF"}{print $0}'> col.4
    rm bed
    
    paste col.* | awk '{printf "%s\t", $1; for(i=2;i<=NF;i+=2) printf "%s\t",\
        $i;print ""}' > $featureDir/$data-4C-Window-Feature-Matrix.txt
    rm col.*
fi

######## Creating Distance Matrix ################
if [ ! -f "$featureDir/$data-Distance.txt" ]
then
    echo -e "Generating distance matrix..."
    cat $train | sort -u > pairs
    awk 'FNR==NR {x[$1];next} ($5 in x)' pairs $ccREs | \
        awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 }' > enhancers
    python $scriptDir/calculate.distance.py $tss enhancers pairs > \
        $featureDir/$data-Distance.txt 
fi


######## Running Random Forest/GBM ################
echo -e "Running Model..."
cat $train | awk '{print $2}' | sort -u  > genes
awk 'FNR==NR {x[$1];next} ($7 in x)' genes $tss | \
awk '{print $1 "\t" $2-500 "\t" $3+500 "\t" $4 "\t" $7 }' > tss

python $scriptDir/gradient.boosting.4C.py $train $featureDir/$data-4C-Enhancer-Feature-Matrix.txt \
    $featureDir/$data-4C-TSS-Feature-Matrix.txt $featureDir/$data-4C-Window-Feature-Matrix.txt \
    $featureDir/$data-Distance.txt tss $data $outputDir $version

rm -r /tmp/moorej3/$SLURM_JOBID-$jid
