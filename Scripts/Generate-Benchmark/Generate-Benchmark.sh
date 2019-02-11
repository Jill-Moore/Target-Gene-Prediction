#!/bin/bash

biosample=$1
linkType=$2
links=$3
name=$4

genome=hg19
version=V4
scriptDir=~/Projects/Target-Gene-Prediction/Scripts/Generate-Benchmark
masterDir=~/Lab/ENCODE/Encyclopedia/$version/Registry/$version-$genome/
masterList=$masterDir/Cell-Type-Specific/Master-Cell-List.txt

if [[ $genome == "mm10" ]]
then
tss=~/Lab/Reference/Mouse/GencodeM4/TSS.Filtered.4K.bed
elif [[ $genome == "hg38" ]]
then
tss=~/Lab/Reference/Human/$genome/GENCODE24/TSS.Filtered.4K.bed
elif [[ $genome == "hg19" ]]
then
tss=~/Lab/Reference/Human/$genome/Gencode19/TSS.Filtered.4K.bed
fi

file=$(awk '{if ($9 == "'$biosample'") print $2}' $masterList)

bedtools intersect -v -a $masterDir/Cell-Type-Specific/Five-Group/$file*.bed \
    -b $tss | grep "Enhancer" > enhancers.bed

if [ $linkType == "ChIA-PET" ]
then
    awk '{if ($NF >= 4) print $0}' $links > links
    python $scriptDir/process.chiapet.py links enhancers.bed $tss \
        $name-Blacklist.txt > $name-Links.txt

elif [ $linkType == "Hi-C" ]
then
    awk '{if (NR != 1) print "chr"$1 "\t" $2 "\t" $3 "\t" "chr"$4 \
        "\t" $5 "\t" $6}' $links > links
    python $scriptDir/process.chiapet.py links enhancers.bed $tss \
        $name-Blacklist.txt > $name-Links.txt
    rm links

elif [ $linkType == "CHi-C" ]
then
    awk '{if (NR != 1 && $NF > 10) print $1 "\t" $2 "\t" $3 "\t" $7 \
        "\t" $8 "\t" $9}' $links > links
    python $scriptDir/process.chiapet.py links enhancers.bed $tss \
        $name-Blacklist.txt > $name-Links.txt
    rm links

elif [ $linkType == "eQTL" ]
then
    python $scriptDir/process.eqtl.py $links enhancers.bed $tss \
        > $name-Links.txt  

elif [ $linkType == "CRISPR" ]
then
    python $scriptDir/process.crispr.py $links enhancers.bed $tss \
        > $name-Links.txt
    tssK562=~/Lab/Target-Gene/Benchmark/Raw/CRISPR/Shendure/K562-V19-TSS.bed
fi

if [ $linkType == "CRISPR" ]
then
    cutoff=1000000
    python $scriptDir/calculate.distance.py $tss enhancers.bed \
        $name-Links.txt $name-Distance.txt
else
    cutoff=$(python $scriptDir/calculate.distance.py $tss enhancers.bed \
        $name-Links.txt $name-Distance.txt)
fi

awk '{if ($3 <= '$cutoff') print $0}' $name-Distance.txt > \
    $name-Distance.cutoff.txt

if [ $linkType == "CRISPR" ]
then
    python $scriptDir/create.experiment.sets.py $name-Distance.cutoff.txt \
        $tssK562 enhancers.bed output $cutoff $name $linkType
else
    python $scriptDir/create.experiment.sets.py $name-Distance.cutoff.txt \
        $tss enhancers.bed output $cutoff $name $linkType
fi

awk '{print $1}' positive | sort -u > ccre-list.txt
awk 'FNR==NR {x[$1];next} ($4 in x)' ccre-list.txt enhancers.bed > tmp.bed
awk '{print $1 "\t" $2 "\t" 1}' positive > total
awk '{print $1 "\t" $2 "\t" 0}' negative >> total

python $scriptDir/assign.groups.py tmp.bed total > $name-Benchmark.v1.txt

rm -f bed1 bed2 out1 out2 positive negative range output
rm -f ccre-list.txt enhancers.bed intersection2.bed total tmp.bed
