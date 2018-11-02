#!/bin/bash

setDir=~/Lab/Target-Gene/Benchmark
scriptDir=~/Projects/Target-Gene-Prediction/Scripts/Benchmark-Characteristics
outputDir=~/Lab/Target-Gene/Benchmark/Characteristics
mkdir -p $outputDir

echo "" > header
datasets=("Aiden-HiC" "Ruan-CTCF" "Ruan-RNAPII" "CHiC" "GEUVADIS" "GTEx")
for data in ${datasets[@]}
do
    echo $data >> header
done
datasets=("Aiden-HiC" "Ruan-CTCF" "Ruan-RNAPII" "CHiC" "GEUVADIS" "GTEx")
for data1 in ${datasets[@]}
do
    echo $data1 > col
    d1=$setDir/$data1-Benchmark.v1.txt 
    for data2 in ${datasets[@]}
    do
        echo -e "\t"$data2
        d2=$setDir/$data2-Benchmark.v1.txt
        python $scriptDir/overlap.coefficient.py $d1 $d2 >> col
    done
    paste header col > tmp
    mv tmp header
done
mv header $outputDir/Benchmark-Overlap-Matrix.txt

rm col
