#!/bin/bash
#___UNDOCUMENTED___

echo "Starting ends analysis run at" `date`

NAME=$1
GENOME=$2
annotations=$3
contexts=$4
orientations=$5
FLAG=$6
OUTDIR=$7

avoid=( groupmt mt chrmt chrc chrm dmel_mitochondrion_genome group_c519-c654 scaffold_m162 )
groups=`grep '>' $GENOME | perl -ne 's/>([^\s]+)//; print $1, "\n"' | tr "[:upper:]" "[:lower:]"`

cd $NAME
#rm -r ${OUTDIR}
echo "Computing scores"
mkdir -p ${OUTDIR}/{log,dat}

for i in $groups; do
    
    skipped=`perl -e '$a = shift; $avoid = join q{|}, @ARGV; if ($a =~ /($avoid)/) {print $1; exit 1} else {exit 0}' $i ${avoid[@]}`
    result=$?

    if [ $result -eq 0 ]; then

        for j in $contexts; do
            for k in $annotations; do
                for l in $orientations; do

                    atype=`echo $k | cut -d'|' -f1`
                    annotation=`echo $k | cut -d'|' -f2 | sed s/\.gff/-$i.gff/`
                    cat_flag=`echo $FLAG | sed 's/ /-/'`
                    fin=${NAME}_BS-Seq_${i}_${j}_w1_methylation.gff
                    fout=${NAME}_BS-Seq_${i}_${j}_${atype}_${l}_flag${cat_flag}
                    makeDist.pl $annotation post-processing/single-c/${fin} ${OUTDIR}/${fout} 100 5000 $FLAG $l
                    mv ${OUTDIR}/${fout}.log ${OUTDIR}/log/${fout}.log
                    average.pl ${OUTDIR}/${fout}.ends 100 100 > ${OUTDIR}/dat/${fout}.dat
                    
                done
            done
        done

    fi
    
done


for j in $contexts; do
    for k in $annotations; do
        for l in $orientations; do
            
            atype=`echo $k | cut -d'|' -f1`
            cat_flag=`echo $FLAG | sed 's/ /-/'`

            cat-dat.pl ${OUTDIR}/dat/${NAME}_BS-Seq_*_${j}_${atype}_${l}_flag${cat_flag}.dat > ${OUTDIR}/dat/${NAME}_BS-Seq_all_per_chr_${j}_${atype}_${l}_flag${cat_flag}.dat
            
        done
    done
done


echo "Finished ends analysis run at" `date`
