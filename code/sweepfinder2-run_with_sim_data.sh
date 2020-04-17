conda activate sch_man_nwinvasion-sweepfinder

cd /master/nplatt/sch_man_nwinvasion/results/sweepfinder
mkdir logs 

CONDA="conda activate sch_man_nwinvasion-sweepfinder"
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 2 "    

for POP in senegal tanzania niger brazil; do
    mkdir $POP-sim

    for PROBED_VCF in $(ls ../msprime/$POP/chr1_"$POP"_rep_*seed_*_probed.vcf); do
        SAMPLE_NAME=$(basename $PROBED_VCF .vcf)
        
        #convert to sweep format
        vcftools \
            --vcf $PROBED_VCF  \
            --counts2 \
            --stdout \
            >$POP-sim/$SAMPLE_NAME.freq

        tail -n+2 $POP-sim/$SAMPLE_NAME.freq  | awk -v OFS="\t" '{print $2,$6,$4,"1"}' >$POP-sim/$SAMPLE_NAME.in

        echo -e 'position\tx\tn\tfolded' | cat - $POP-sim/$SAMPLE_NAME.in > temp && mv temp $POP-sim/$SAMPLE_NAME.in

        #generate sfs
        SFS_CMD="SweepFinder2 -f $POP-sim/$SAMPLE_NAME.in $POP-sim/$SAMPLE_NAME.sfs >logs/"$SAMPLE_NAME"_sfs.log"

        #gen sweepfinder2
        SF2_CMD="SweepFinder2 -lg 1000 $POP-sim/$SAMPLE_NAME.in $POP-sim/$SAMPLE_NAME.sfs $POP-sim/$SAMPLE_NAME.sw2out"
        
        #submit to the queue
        CMD="$CONDA; $SFS_CMD; $SF2_CMD"
        echo $CMD | $QSUB -N sf2s_$SAMPLE_NAME -o logs/sf2_sim_$SAMPLE_NAME.log

    done
done

