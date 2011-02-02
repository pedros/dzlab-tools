#!/bin/bash
#___UNDOCUMENTED___

NAME=$1
GENOME=$2
groups=`grep '>' $GENOME | perl -ne 's/>([^\s]+)//; print $1, "\n"' | tr "[:upper:]" "[:lower:]" `
avoid=('groupmt' 'mt' 'chrmt' 'chrc' 'chrm' 'dmel_mitochondrion_genome' 'group_c519-c654' 'scaffold_m162')

echo "Working on $NAME"

cd $NAME
rm    -r post-processing/freq
mkdir -p post-processing/freq

for i in b*; do
    for j in ${i}/single-c/*freq; do

        skipped=`perl -e '$a = shift; $avoid = join q{|}, @ARGV; if ($a =~ /($avoid)/) {print $1; exit 1} else {exit 0}' $j ${avoid[@]}`
        result=$?
        
        if [ $result -eq 0 ]; then
            cat $j >> post-processing/freq/${NAME}_BS-Seq_nuclear_${i}_methylation.freq
        else
            cat $j >> post-processing/freq/${NAME}_BS-Seq_${skipped}_${i}_methylation.freq
        fi
    done
    merge_freq.pl post-processing/freq/${NAME}_BS-Seq_nuclear_${i}_methylation.freq >> post-processing/freq/${NAME}_BS-Seq_nuclear_per-batch_methylation.freq
    rm post-processing/freq/${NAME}_BS-Seq_nuclear_${i}_methylation.freq
done

for i in $groups; do
    skipped=`perl -e '$a = shift; $avoid = join q{|}, @ARGV; if ($a =~ /($avoid)/) {print $1; exit 1} else {exit 0}' $i ${avoid[@]}`
    result=$?

    for j in b*; do
        if [ $result -eq 0 ]; then
            cat ${j}/single-c/*[-_.]${i}[-_.]*freq >> post-processing/freq/${NAME}_BS-Seq_nuclear_methylation.freq
        else
            cat ${j}/single-c/*[-_.]${i}[-_.]*freq >> post-processing/freq/${NAME}_BS-Seq_${skipped}_methylation.freq
            cat post-processing/freq/${NAME}_BS-Seq_${skipped}_${j}_methylation.freq >> post-processing/freq/${NAME}_BS-Seq_${skipped}_per-batch_methylation.freq
            rm post-processing/freq/${NAME}_BS-Seq_${skipped}_${j}_methylation.freq
        fi
    done
    
    if [ $result -eq 1 ]; then
        merge_freq.pl post-processing/freq/${NAME}_BS-Seq_${skipped}_methylation.freq > post-processing/freq/${NAME}_BS-Seq_${skipped}_methylation.freq.new
        mv post-processing/freq/${NAME}_BS-Seq_${skipped}_methylation.freq.new post-processing/freq/${NAME}_BS-Seq_${skipped}_methylation.freq
    fi
done

merge_freq.pl post-processing/freq/${NAME}_BS-Seq_nuclear_methylation.freq > post-processing/freq/${NAME}_BS-Seq_nuclear_methylation.freq.new
mv post-processing/freq/${NAME}_BS-Seq_nuclear_methylation.freq.new post-processing/freq/${NAME}_BS-Seq_nuclear_methylation.freq

parse_fasta.pl -l $GENOME > post-processing/freq/bp_counts

for i in $groups; do
    skipped=`perl -e '$a = shift; $avoid = join q{|}, @ARGV; if ($a =~ /($avoid)/) {print $1; exit 1} else {exit 0}' $i ${avoid[@]}`
    result=$?
    if [ $result -eq 1 ]; then
        let non_nuclear_bp=`cut -f1 '	' post-processing/freq/${NAME}_BS-Seq_${skipped}_methylation.freq | tail -n1 `
        let genome_size=`grep $skipped post-processing/freq/bp_counts  | cut -f3 -d'	'`
        coverage="`perl -e 'printf("%.1f", $ARGV[0]/$ARGV[1])' $non_nuclear_bp $genome_size`"
        echo "Genome size:	$genome_size" >> post-processing/freq/${skipped}_coverage
        echo "Aligned bps:	$non_nuclear_bp" >> post-processing/freq/${skipped}_coverage
        echo "Coverage:	$coverage" >> post-processing/freq/${skipped}_coverage
    fi
done    

let nuclear_bp=`cut -f1 '	' post-processing/freq/${NAME}_BS-Seq_nuclear_methylation.freq | tail -n1`
let genome_size=`grep 'Total size' post-processing/freq/bp_counts  | cut -f3 -d'	'`
coverage="`perl -e 'printf("%.1f", $ARGV[0]/$ARGV[1])' $nuclear_bp $genome_size`"

echo "Genome size:	$genome_size" >> post-processing/freq/nuclear_coverage
echo "Aligned bps:	$nuclear_bp" >> post-processing/freq/nuclear_coverage
echo "Coverage:	$coverage" >> post-processing/freq/nuclear_coverage

cd -
