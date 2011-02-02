#!/bin/sh
#___UNDOCUMENTED___

echo "Starting cat and window run at" `date`

NAME=$1
GENOME=$2
contexts=$3
shift
shift
shift
BATCHES=$@
groups=`grep '>' $GENOME | perl -ne 's/>([^\s]+)//; print $1, "\n"' | tr "[:upper:]" "[:lower:]" `
width=50
step=50

cd $NAME
# rm -r post-processing/single-c
echo "Concatenating single c files..."
mkdir -p post-processing/single-c
for i in $groups; do
    for j in $BATCHES; do
        for k in $contexts; do
            cat ${j}/single-c/*[._-]${i}[._-]single-c*${k}*gff >> post-processing/single-c/${NAME}_BS-Seq_${i}_${k}_w1_methylation.gff
        done
    done
    for k in $contexts; do
        window_gff.pl -f post-processing/single-c/${NAME}_BS-Seq_${i}_${k}_w1_methylation.gff -w 1 -s 1 -o post-processing/single-c/${NAME}_BS-Seq_${i}_${k}_w1_methylation.gff.merged
        mv post-processing/single-c/${NAME}_BS-Seq_${i}_${k}_w1_methylation.gff.merged post-processing/single-c/${NAME}_BS-Seq_${i}_${k}_w1_methylation.gff
    done
done
echo "Done with code: $?"


echo "Merging and windowing concatenated single c files"
# rm -r post-processing/windows
mkdir -p post-processing/windows
for i in $groups; do
    for j in $contexts; do
        window_gff_REFACTORED.pl post-processing/single-c/${NAME}_BS-Seq_${i}_${j}_w1_methylation.gff --width ${width} --step ${step} --absolute rice --no-skip --output post-processing/windows/${NAME}_BS-Seq_${i}_${j}_w${width}_methylation.gff
        cat post-processing/windows/${NAME}_BS-Seq_${i}_${j}_w${width}_methylation.gff >> post-processing/windows/${NAME}_BS-Seq_all_${j}_w${width}_methylation.gff
    done
done
echo "Done with code: $?"
cd ..



echo "Finished cat and window run at" `date`
