#!/bin/sh

echo "Starting cat and window run at" `date`

NAME=$1
GENOME=$2
shift
shift
BATCHES=$@
groups=`grep '>' $GENOME | perl -ne 's/>([^\s]+)//; print $1, "\n"' | tr "[:upper:]" "[:lower:]" `
width=50
step=50
contexts='CG CHG CHH'

# GSM399598_WT_endosperm_BS_seq_w50-wt_endosperm_CHH

echo "Concatenating single c files..."
mkdir -p post-processing/single-c
for i in $groups; do
    for j in $BATCHES; do
        for k in $contexts; do
            # ls ${j}/single-c/*${i}*single-c*${k}*gff
            # ls post-processing/single-c/${NAME}_BS-Seq_${i}_${k}_w1_methylation.gff
            cat ${j}/single-c/*${i}*single-c*${k}*gff >> post-processing/single-c/${NAME}_BS-Seq_${i}_${k}_w1_methylation.gff
        done
    done
done
echo "Done with code: $?"


echo "Merging and windowing concatenated single c files"
mkdir -p post-processing/windows
for i in $groups; do
    for j in $contexts; do

        # ls post-processing/single-c/${NAME}_BS-Seq_${i}_${k}_w1_methylation.gff
        # ls post-processing/windows/${NAME}_BS-Seq_${i}_${k}_w${width}_methylation.gff
        window_gff.pl --gff-file post-processing/single-c/${NAME}_BS-Seq_${i}_${k}_w1_methylation.gff --width ${width} --step ${step} --no-skip --output post-processing/windows/${NAME}_BS-Seq_${i}_${k}_w${width}_methylation.gff

    done
done
echo "Done with code: $?"


echo "Finished cat and window run at" `date`
