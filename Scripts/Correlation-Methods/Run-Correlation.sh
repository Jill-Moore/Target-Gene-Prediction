


signal=H3K27ac
mode=raw
scriptDir=~/Projects/Target-Gene-Prediction/Scripts/Correlation-Methods
featureDir=~/Lab/Target-Gene/Correlation-Methods/Signal-Matrices

######## Creating Enhancer Signal Matrix ################
if [ ! -f "$featureDir/ccRE-$signal-Matrix.$mode.txt" ]
then
    echo -e "Generating ccRE signal matrix..."
    fileDir=~/Lab/ENCODE/Encyclopedia/V4/Registry/V4-hg19/
    file=$fileDir/$signal-List.txt
    mkdir -p $fileDir/signal-output/$signal
    q=$(wc -l $file | awk '{print $1}')
    for j in `seq 1 1 $q`
    do
        echo $j
        A=$(awk '{if (NR == '$j') print $1}' $file)
        B=$(awk '{if (NR == '$j') print $2}' $file)
        C=$( awk '{if (NR == '$j') print $3}' $file)
        mv $fileDir/signal-output/$A"-"$B.txt $fileDir/signal-output/$signal/
    done
    cd $fileDir/signal-output/$signal/
    if [ $mode == "raw" ]
    then
        paste *.txt | awk '{printf "%s\t", $1; for(i=3;i<=NF;i+=4) printf \
            "%s\t",$i ; print ""}' > matrix
    else
        paste *.txt | awk '{printf "%s\t", $1; for(i=2;i<=NF;i+=4) printf \
            "%s\t",$i ; print ""}' > matrix
    fi
    cp matrix $featureDir/ccRE-$signal-Matrix.$mode.txt
    mv *.txt ../
fi


if [ ! -f "$featureDir/TSS-$signal-Matrix.$mode.txt" ]
then
    echo -e "Generating TSS signal matrix..."
    fileDir=~/Lab/ENCODE/Encyclopedia/V4/Registry/V4-hg19/
    file=$fileDir/$signal-List.txt
    output=$featureDir/signal-output
    tss=~/Lab/Reference/Human/hg19/Gencode19/TSS.Filtered.bed 
    
    mkdir -p $output
    
    if [ $signal == "DNase" ]
    then
        width=250
    elif [ $signal == "H3K27ac" ]
    then
        width=750
    fi
    num=$(wc -l $file | awk '{print $1}')
    jobid=$(sbatch --nodes 1 --array=1-$num%50 --mem=5G --time=04:00:00 \
        --output=/home/moorej3/Job-Logs/jobid_%A_%a.output \
        --error=/home/moorej3/Job-Logs/jobid_%A_%a.error \
        $scriptDir/Retrieve-Signal.sh $tss $signal $file $output $width | \
        awk '{print $4}')
    echo $jobid
    sleep 20
    list=100
    while [ $list -gt 1 ]
    do
        list=$(squeue -j $jobid | wc -l | awk '{print $1}')
        echo -e "jobs still running: $list"
        sleep 10
    done
    
    cd $output
    if [ $mode == "raw" ]
    then
        paste $signal.*.txt | awk '{printf "%s\t", $1; for(i=3;i<=NF;i+=4) printf \
            "%s\t",$i ; print ""}' > matrix
    else
        paste $signal.*.txt | awk '{printf "%s\t", $1; for(i=2;i<=NF;i+=4) printf \
            "%s\t",$i ; print ""}' > matrix
    fi
    cp matrix $featureDir/TSS-$signal-Matrix.$mode.txt
fi  
    

