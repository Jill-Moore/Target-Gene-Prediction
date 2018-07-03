#!/bin/bash
#Jill E. Moore
#Weng Lab
#UMass Medical School
#Updated May 2018

#CTCF, DNase, H3K27ac, H3K4me3
mode=$2
genome=$1

fileDir=~/Lab/ENCODE/Encyclopedia/V4/Registry/V4-$genome

echo "Processing Files ..."
file=$fileDir/$mode-List.txt
mkdir -p $fileDir/signal-output/$mode
q=$(wc -l $file | awk '{print $1}')
for j in `seq 1 1 $q`
do
    A=$(awk '{if (NR == '$j') print $1}' $file)
    B=$(awk '{if (NR == '$j') print $2}' $file)
    C=$( awk '{if (NR == '$j') print $3}' $file)
    mv $fileDir/signal-output/$A"-"$B.txt $fileDir/signal-output/$mode/
done

echo "Creating Matrix..."
cd $fileDir/signal-output/$mode/
paste *.txt | \
awk '{printf "%s\t", $1; for(i=2;i<=NF;i+=4) printf "%s\t",$i ; print ""}' > matrix

echo "Determining maxZ..."
python ~/Projects/Target-Gene-Prediction/Scripts/count-biosamples.py  matrix > \
    $fileDir/ccRE-$mode-Biosample-Counts.txt
mv *.txt ../

