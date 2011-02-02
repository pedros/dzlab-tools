#!/bin/bash
#___UNDOCUMENTED___

echo "Starting run at " `date`

LREAD=$1
RREAD=$2
GENOME=$3
TISSUE=$4
BATCH=$5
OVERWRITE=$6

if [ ! -e ${LREAD}.fa -o $OVERWRITE == 1 ]; then
    # convert from fastq to fasta
    echo -n "Converting sequences to fasta format..."
    fq_all2std.pl fq2fa $LREAD > ${LREAD}.fa
    echo "done with code: $?"
fi

if [ ! -e ${RREAD}.fa -o $OVERWRITE == 1 ]; then
    # convert from fastq to fasta
    echo -n "Converting sequences to fasta format..."
    fq_all2std.pl fq2fa $RREAD > ${RREAD}.fa
    echo "done with code: $?"
fi

if [ ! -e ${LREAD}.c2t -o $OVERWRITE == 1 ]; then
    # convert sequences
    echo -n "Converting sequences..."
    convert.pl c2t ${LREAD}.fa > ${LREAD}.c2t
    echo "done with code: $?"
fi

if [ ! -e ${RREAD}.c2t -o $OVERWRITE == 1 ]; then
    # convert sequences
    echo -n "Converting sequences..."
    convert.pl g2a ${RREAD}.fa > ${RREAD}.g2a
    echo "done with code: $?"
fi


if [ ! -e ${GENOME}-RC.G2A -o $OVERWRITE == 1 ]; then
    # convert genome
    echo -n "Converting genome..."
    rcfas.pl $GENOME > ${GENOME}-RC
    convert.pl c2t ${GENOME}-RC > ${GENOME}-RC.C2T
    convert.pl g2a ${GENOME}-RC > ${GENOME}-RC.G2A
    echo "done with code: $?"
fi

# /available_memory:`free -m | grep 'Mem:' | perl -e 'while(<>) {@a=(split /\s+/, $_); print $a[3];}'`
# align sequences
echo -n "Aligning sequences..."
seqmap 2 ${LREAD}.c2t ${GENOME}-RC.C2T ${LREAD}.eland3 /eland:3 /forward_strand /available_memory:8000 /cut:1,45 &&
seqmap 2 ${RREAD}.g2a ${GENOME}-RC.G2A ${RREAD}.eland3 /eland:3 /forward_strand /available_memory:8000 /cut:1,45 &&
# echo "done with code: $?"

# correlate paired ends
echo -n "Running correlatePairedEnds.pl..."
correlatePairedEnds.pl --left ${LREAD}.eland3 --right ${RREAD}.eland3 --reference $GENOME --output ${LREAD}_pre.gff --offset 0 --distance 300 --readsize 45 &&
echo "Done with code: $?"

# replace processed reads with original
echo -n "Running replaceMutation.pl..."
replaceMutation.pl ${LREAD}.fa ${RREAD}.fa ${LREAD}_pre.gff 45 > ${LREAD}_post.gff &&
echo "Done with code: $?"

# check alignment ${LREAD}.accuracy
echo -n "Computing alignment ${LREAD}.accuracy..."
collect_align_stats.pl ${LREAD}.eland3 ${RREAD}.eland3 ${LREAD}_pre.gff $TISSUE $BATCH > ${LREAD}_pre.alignment.log
echo "Done with code: $?"

# filter out all non matches
echo -n "Filtering out records without matches..."
grep 'target=' ${LREAD}_post.gff > ${LREAD}_post_filtered.gff &&
echo "Done with code: $?"

# split into multiple chromosomes
echo -n "Splitting into multiple chromosomes..."
for i in `grep '>' $GENOME | sed s/\>//`
do
    grep -i "^$i	" ${LREAD}_post_filtered.gff > ${LREAD}_post_filtered_${i}.gff
done
echo "Done with code: $?"

# count methylation
echo -n "Running countMethylation.pl..."
for i in `grep '>' $GENOME | sed s/\>//`
do
    countMethylation.pl --ref $GENOME --gff ${LREAD}_post_filtered_${i}.gff --output ${LREAD}_post_filtered_${i}_singleC.gff --sort
done
echo "Done with code: $?"

# split into multiple contexts
echo -n "Splitting into multiple contexts..."
for i in `grep '>' $GENOME | sed s/\>//`
do
    for j in CG CHG CHH
    do
        grep $j ${LREAD}_post_filtered_${i}_singleC.gff > ${LREAD}_post_filtered_${i}_${j}_singleC.gff
    done
done
echo "Done with code: $?"

# # window files
# echo "Windowing files..."
# window_gff.pl -b ${LREAD}_*_singleC.gff -w 50 -s 50
# echo "Done with code: $?"


## Done
echo "Finished run at " `date`
