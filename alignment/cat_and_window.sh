#!/bin/sh

echo "Starting cat and window run at" `date`

NAME=$1
GENOME=$2
GENES=$3
TRANSPOSONS=$4
contexts=$5
shift
shift
shift
shift
shift
BATCHES=$@
groups=`grep '>' $GENOME | perl -ne 's/>([^\s]+)//; print $1, "\n"' | tr "[:upper:]" "[:lower:]" `
width=50
step=50

echo "Concatenating single c files..."
mkdir -p post-processing/single-c
for i in $groups; do
    for j in $BATCHES; do
        for k in $contexts; do
            # ls ${j}/single-c/*${i}*single-c*${k}*gff
            # ls post-processing/single-c/${NAME}_BS-Seq_${i}_${k}_w1_methylation.gff
            cat ${j}/single-c/*${i}*single-c*${k}*gff >> post-processing/single-c/${NAME}_BS-Seq_${i}_${k}_w1_methylation.gff
            cat ${j}/single-c/*${i}*single-c.gff.freq >> post-processing/single-c/${NAME}_BS-Seq_${i}_methylation.freq
        done
    done
done
echo "Done with code: $?"

# echo "split_gff.pl  --sequence all $GENES"
# echo "split_gff.pl  --sequence all $TRANSPOSONS"

echo "Merging and windowing concatenated single c files"
mkdir -p post-processing/windows
mkdir -p post-processing/ends-analysis/dat
for i in $groups; do
    for j in $contexts; do

        # ls post-processing/single-c/${NAME}_BS-Seq_${i}_${k}_w1_methylation.gff
        # ls post-processing/windows/${NAME}_BS-Seq_${i}_${k}_w${width}_methylation.gff
        window_gff.pl --gff-file post-processing/single-c/${NAME}_BS-Seq_${i}_${j}_w1_methylation.gff --width ${width} --step ${step} --no-skip --output post-processing/windows/${NAME}_BS-Seq_${i}_${j}_w${width}_methylation.gff
        
        # for annotation in genes transposons; do
        #     for origin in 3prime 5prime; do
        #         for wflag in 2 '6-1500'; do

        #             if [ "$annotation" = "genes" ]; then
        #                 ANNOTATION=$GENES
        #             else
        #                 ANNOTATION=$TRANSPOSONS
        #             fi
        #             flag=`echo $wflag | sed 's/-/ /'`

        #             echo "makeDist.pl `echo $ANNOTATION | sed 's/\.gff/-${i}.gff/'` post-processing/single-c/${NAME}_BS-Seq_${i}_${j}_w1_methylation.gff post-processing/ends-analysis/${NAME}_BS-Seq_${i}_${j}_w1_methylation_${annotation}_${origin}_${wflag} 100 5000 $flag $origin"
        #             echo "average.pl post-processing/ends-analysis/${NAME}_BS-Seq_${i}_${j}_w1_methylation_${annotation}_${origin}_${flag} 100 100 > post-processing/ends-analysis/dat/${NAME}_BS-Seq_${i}_${j}_w1_methylation_${annotation}_${origin}_${wflag}.dat"
        #         done
        #     done
        # done

    done
done
echo "Done with code: $?"




echo "Finished cat and window run at" `date`
