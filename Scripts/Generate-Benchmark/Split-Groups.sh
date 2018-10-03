data=$1

awk '{print $1}' $data-Positive.txt | sort -u > cre-list.txt
half=$(wc -l cre-list.txt | awk '{printf "%.0f", $1/2}')

sort -R cre-list.txt | head -n $half > train
awk 'FNR==NR {x[$0];next} !($0 in x)' train cre-list.txt > tmp
half=$(wc -l tmp | awk '{printf "%.0f", $1/2}')

sort -R tmp | head -n $half > val
awk 'FNR==NR {x[$0];next} !($0 in x)' val tmp > test

awk 'FNR==NR {x[$1];next} ($1 in x)' train $data-Positive.txt > $data-Training-P.txt
awk 'FNR==NR {x[$1];next} ($1 in x)' train $data-Negative.txt > $data-Training-N.txt

awk 'FNR==NR {x[$1];next} ($1 in x)' test $data-Positive.txt > $data-Test-P.txt
awk 'FNR==NR {x[$1];next} ($1 in x)' test $data-Negative.txt > $data-Test-N.txt

awk 'FNR==NR {x[$1];next} ($1 in x)' val $data-Positive.txt > $data-Validation-P.txt
awk 'FNR==NR {x[$1];next} ($1 in x)' val $data-Negative.txt > $data-Validation-N.txt

