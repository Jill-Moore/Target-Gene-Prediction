#!/bin/bash

data=$1
biosample=$(echo $data | awk -F "." '{print $1}')

setDir=~/Lab/Target-Gene/Benchmark
scriptDir=~/Projects/Target-Gene-Prediction/Scripts/Benchmark-Characteristics
train=$setDir/$data-Benchmark.v3.txt
output=~/Lab/Target-Gene/Benchmark/Characteristics/TF-Signal
dataDir=/data/projects/encode/data/
tss=~/Lab/Reference/Human/hg19/Gencode19/TSS.2019.bed
ccres=~/Lab/ENCODE/Encyclopedia/V4/Registry/V4-hg19/hg19-ccREs-Simple.bed
tfList=~/Lab/Target-Gene/Benchmark/Characteristics/CTCF-List.txt

cat $train | awk '{print $1}' | sort -u  > ccres
awk 'FNR==NR {x[$1];next} ($5 in x)' ccres $ccres > tmp1

#cat $train | awk '{print $2}' | sort -u  > genes
#awk 'FNR==NR {x[$1];next} ($7 in x)' genes $tss | \
#    awk '{print $1 "\t" $2-250 "\t" $3+250 "\t" $4}' > tss

width=0

dset=$(grep $biosample $tfList | awk '{print $1}')
dsig=$(grep $biosample $tfList | awk '{print $2}')

awk -F "\t" '{printf "%s\t%.0f\t%.0f\t%s\n", $1,$2-'$width',$3+'$width',$4}' \
 tmp1 | awk '{if ($2 < 0) print $1 "\t" 0 "\t" $3 "\t" $4 ; else print $0}' \
| sort -u > little

~/bin/bigWigAverageOverBed -bedOut=out2.bed $dataDir/$dset/$dsig.bigWig little out2

awk '{print $1 "\t" $5}' out2 | sort -k1,1 > $output/$data"-ELS."$dset"-"$dsig".txt" 

rm ccres tmp1 little out2.bed out2
