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

grep "Enhancer" $masterDir/Cell-Type-Specific/Five-Group/$file*.bed > \
    enhancers.bed


if [ $linkType == "ChIA-PET" ] || [ $linkType == "Hi-C" ]
then
    python $scriptDir/process.chiapet.py $links enhancers.bed $tss \
        $name-Blacklist.txt > $name-Links.txt
elif [ $linkType == "eQTL" ]
then
    python $scriptDir/process.chiapet.py $links enhancers.bed $tss \
        $name-Blacklist.txt > $name-Links.txt  
fi

cutoff=$(python $scriptDir/calculate.distance.py $tss enhancers.bed \
    $name-Links.txt $name-Distance.txt)

awk '{if ($3 >= '$cutoff') print $0}' $name-Distance.txt > \
    $name-Distance.cutoff.txt

python $scriptDir/create.experiment.sets.py $name-Distance.cutoff.txt \
    $tss enhancers.bed output $cutoff $name-Blacklist.txt

awk '{print $1}' positive | sort -u > ccre-list.txt
half=$(wc -l ccre-list.txt | awk '{printf "%.0f", $1/2}')

sort -R ccre-list.txt | head -n $half > train
awk 'FNR==NR {x[$0];next} !($0 in x)' train ccre-list.txt > tmp
half=$(wc -l tmp | awk '{printf "%.0f", $1/2}')

sort -R tmp | head -n $half > val
awk 'FNR==NR {x[$0];next} !($0 in x)' val tmp > test

awk 'FNR==NR {x[$1];next} ($1 in x)' train positive | \
    awk '{print $1 "\t" $2 "\t" 1}' > $name-Training.txt
awk 'FNR==NR {x[$1];next} ($1 in x)' train negative | \
    awk '{print $1 "\t" $2 "\t" 0}' >> $name-Training.txt

awk 'FNR==NR {x[$1];next} ($1 in x)' test positive | \
    awk '{print $1 "\t" $2 "\t" 1}' > $name-Test.txt
awk 'FNR==NR {x[$1];next} ($1 in x)' test negative | \
    awk '{print $1 "\t" $2 "\t" 0}' >> $name-Test.txt

awk 'FNR==NR {x[$1];next} ($1 in x)' val positive | \
    awk '{print $1 "\t" $2 "\t" 1}' > $name-Validation.txt
awk 'FNR==NR {x[$1];next} ($1 in x)' val negative | \
    awk '{print $1 "\t" $2 "\t" 0}' >> $name-Validation.txt

cat $name-Training.txt $name-Test.txt $name-Validation.txt > $name-Total.txt

rm bed1 bed2 out1 out2 tmp train val positive negative range output test
rm ccre-list.txt enhancers.bed intersection2.bed
