
#!/bin/bash
#SBATCH --nodes 1
#SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH --output=/home/moorej3/Job-Logs/jobid_%A.output
#SBATCH --error=/home/moorej3/Job-Logs/jobid_%A.error
#SBATCH --partition=12hours

data=Aiden-HiC
version=v1

setDir=~/Lab/Target-Gene/Benchmark
train=$setDir/$data-Benchmark.$version.txt

sigDir=/data/projects/roadmap/data/imputed/

featureDir=~/Lab/Target-Gene/JEME/Feature-Matrices
outputDir=~/Lab/Target-Gene/JEME/Results
enhancers=~/Lab/Target-Gene/Target-Finder/GM12878-Enhancers.bed
scriptDir=~/Projects/Target-Gene-Prediction/Scripts/JEME
tss=~/Lab/Reference/Human/hg19/Gencode19/TSS.Filtered.bed
regEl=~/Lab/ENCODE/Encyclopedia/V4/Registry/V4-hg19/hg19-ccREs-Simple.bed

mkdir -p /tmp/moorej3/$SLURM_JOBID-$jid
cd /tmp/moorej3/$SLURM_JOBID-$jid

mkdir -p $outputDir

######## Creating Enhancer Feature Matrix ################
if [ ! -f "$featureDir/$data-Enhancer-Feature-Matrix.txt" ]
then
    echo -e "Generating enhancer feature matrix..."
    cat $train | awk '{print $1}' | sort -u  > cres
    for sig in DNase H3K4me1 H3K27ac H3K27me3
    do
        bw=$sigDir/$sig/E116-$sig.imputed.pval.signal.bigwig
        if [ $sig == "DNase" ]
        then
            width=0
        elif [ $sig == "H3K27ac" ] || [ $sig == "RNAseq" ] || \
            [ $sig == "H3K4me1" ] || [ $sig == "H3K27me3" ]
        then
            width=500
        fi
        awk 'FNR==NR {x[$1];next} ($4 in x)' cres $enhancers | awk '{print $1 \
            "\t" $2-'$width' "\t" $3+'$width' "\t" $4}' > enhancers
        head -n 1 enhancers
        ~/bin/bigWigAverageOverBed $bw enhancers out
        sort -k1,1 out | awk 'BEGIN{print "ccre\t'$sig'"}{print $1 "\t" $5}' \
            > out.$sig
    done
    paste out.DNase out.H3K4me1 out.H3K27ac out.H3K27me3 | awk  '{printf \
        "%s", $1; for(i=2;i<=NF;i+=2) printf "\t%s",$i ; print ""}'  > \
        $featureDir/$data-Enhancer-Feature-Matrix.txt
    rm out.* enhancers
fi

######## Creating TSS Feature Matrix ################
if [ ! -f "$featureDir/$data-TSS-Feature-Matrix.txt" ]
then
    echo -e "Generating tss feature matrix..."
    cat $train | awk '{print $2}' | sort -u  > genes
    awk 'FNR==NR {x[$1];next} ($7 in x)' genes $tss | \
        awk '{print $1 "\t" $2-500 "\t" $3+500 "\t" $4 "\t" $7 }' > tss
    for sig in DNase H3K4me1 H3K27ac H3K27me3
    do
        bw=$sigDir/$sig/E116-$sig.imputed.pval.signal.bigwig
        if [ $sig == "DNase" ]
        then
            width=250
        elif [ $sig == "H3K27ac" ] || [ $sig == "RNAseq" ] || \
            [ $sig == "H3K4me1" ] || [ $sig == "H3K27me3" ]
        then
            width=750
        fi
        awk 'FNR==NR {x[$1];next} ($7 in x)' genes $tss | awk '{print $1 \
            "\t" $2-'$width' "\t" $3+'$width' "\t" $4}' > tss
        head -n 1 tss
        ~/bin/bigWigAverageOverBed $bw tss out
        sort -k1,1 out | awk 'BEGIN{print "tss\t'$sig'"}{print $1 "\t" $5}' \
            > out.$sig
    done
    paste out.DNase out.H3K4me1 out.H3K27ac out.H3K27me3 | awk  '{printf \
        "%s", $1; for(i=2;i<=NF;i+=2) printf "\t%s",$i ; print ""}'  > \
        $featureDir/$data-TSS-Feature-Matrix.txt
    rm out.* tss
fi

######## Creating Window Matrix ################
if [ ! -f "$featureDir/$data-Window-Feature-Matrix.txt" ]
then
    echo -e "Generating window feature matrix..."
    cat $train | sort -u > pairs
    python $scriptDir/create.window.py $tss $enhancers pairs | awk '{print $1 \
        "\t" $2 "\t" $3 "\t" $4}'> windows
    
    for sig in DNase H3K4me1 H3K27ac H3K27me3
    do
        bw=$sigDir/$sig/E116-$sig.imputed.pval.signal.bigwig
        if [ $sig == "DNase" ]
        then
            width=250
        elif [ $sig == "H3K27ac" ] || [ $sig == "RNAseq" ] || \
            [ $sig == "H3K4me1" ] || [ $sig == "H3K27me3" ]
        then
            width=750
        fi
        head -n 1 windows
        ~/bin/bigWigAverageOverBed $bw windows out
        sort -k1,1 out | awk 'BEGIN{print "tss\t'$sig'"}{print $1 "\t" $5}' \
            > out.$sig
    done
    paste out.DNase out.H3K4me1 out.H3K27ac out.H3K27me3 | awk  '{printf \
        "%s", $1; for(i=2;i<=NF;i+=2) printf "\t%s",$i ; print ""}'  > \
        $featureDir/$data-Window-Feature-Matrix.txt
    
    
######## Creating Distance Matrix ################
if [ ! -f "$featureDir/$data-Distance.txt" ]
then
    echo -e "Generating distance matrix..."
    cat $train | sort -u > pairs
    python $scriptDir/calculate.distance.py $tss $enhancers pairs > \
        $featureDir/$data-Distance.txt 
fi

######## Creating Error Matrix ################
if [ ! -f "$featureDir/$data-Window-Feature-Matrix.txt" ]
then
    echo -e "Generating error feature matrix..."
    cat $train | sort -u > pairs
    rm -f *.running.out
    
    awk 'BEGIN{printf "%s-%s", "DNase",1; for(i=2;i<=127;i+=1) printf "\t%s-%s",\
        "DNase",i;print ""}' > DNase.running.out
    awk 'BEGIN{printf "%s-%s", "H3K4me1",1; for(i=2;i<=127;i+=1) printf "\t%s-%s",\
        "H3K4me1",i;print ""}' > H3K4me1.running.out
    awk 'BEGIN{printf "%s-%s", "H3K27ac",1; for(i=2;i<=127;i+=1) printf "\t%s-%s",\
        "H3K27ac",i;print ""}' > H3K27ac.running.out
    awk 'BEGIN{printf "%s-%s", "H3K27me3",1; for(i=2;i<=127;i+=1) printf "\t%s-%s",\
        "H3K27me3",i;print ""}' > H3K27me3.running.out
    
    awk 'BEGIN{print "pair"}' > header
    i=$(wc -l pairs | awk '{print $1}')
    for j in `seq 1 1 $i`
    do
        echo $j
        ccre=$(awk '{if (NR == '$j') print $1}' pairs)
        rdhs=$(grep $ccre $regEl | awk '{print $4}')
        gene=$(awk '{if (NR == '$j') print $2}' pairs)
        grep $gene $tss > tss
        x=$(wc -l tss | awk '{print $1}')
        for y in `seq 1 1 $x`
        do
            t=$(awk '{if (NR == '$j') print $4}' tss)
            awk 'BEGIN{print "'$rdhs'-'$t'"}' >> header
            for sig in DNase H3K4me1 H3K27ac H3K27me3
            do
                grep $rdhs $featureDir/Error-Terms/$sig.$t.csv | awk -F "," \
                    '{printf "%s", $2; for(i=3;i<=NF;i+=1) printf "\t%s",$i ; \
                    print ""}' >> $sig.running.out
            done
        done
    done
    paste header *.running.out > $featureDir/$data-Error-Feature-Matrix.txt
fi


######## Running Random Forest/GBM ################
#echo -e "Running Model..."
#cat $train | awk '{print $2}' | sort -u  > genes
#awk 'FNR==NR {x[$1];next} ($7 in x)' genes $tss | \
#awk '{print $1 "\t" $2-500 "\t" $3+500 "\t" $4 "\t" $7 }' > tss

#python $scriptDir/gradient.boosting.py $train $featureDir/$data-Enhancer-Feature-Matrix.txt \
#    $featureDir/$data-TSS-Feature-Matrix.txt $featureDir/$data-Window-Feature-Matrix.txt \
#    tss $data $outputDir $version

rm -r /tmp/moorej3/$SLURM_JOBID-$jid
