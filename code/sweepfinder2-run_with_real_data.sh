conda activate sch_man_nwinvasion-sweepfinder

cd /master/nplatt/sch_man_nwinvasion/results/sweepfinder


for POP in new_world east_africa west_africa; do
    mkdir $POP

    #get the pop specific sfs
    vcftools \
        --keep ~/sch_man_nwinvasion/results/lists/$POP.list \
        --vcf ../variant_filtration/smv7_ex_autosomes.vcf  \
        --counts2 \
        --stdout \
        >$POP/$POP.freq

    tail -n+2 $POP/$POP.freq  | awk -v OFS="\t" '{print $2,$6,$4,"1"}' >$POP/$POP.in

    echo -e 'position\tx\tn\tfolded' | cat - $POP/$POP.in > temp && mv temp $POP/$POP.in

    SweepFinder2 -f $POP/$POP.in $POP/$POP.sfs
 
    #run sweepfinder for each chrom
    for I in $(seq 1 7); do
        CHR=SM_V7_$I

        #convert vcf to sw input format
        vcftools \
             --keep ~/sch_man_nwinvasion/results/lists/$POP.list \
            --vcf ../variant_filtration/smv7_ex_autosomes.vcf \
            --chr $CHR \
            --counts2 \
            --stdout \
            >$POP/$CHR"_"$POP.freq

        tail -n+2 $POP/$CHR"_"$POP.freq  | awk -v OFS="\t" '{print $2,$6,$4,"1"}' >$POP/$CHR"_"$POP.in

        echo -e 'position\tx\tn\tfolded' | cat - $POP/$CHR"_"$POP.in > temp && mv temp $POP/$CHR"_"$POP.in

        #submit sweepfinder        
        CMD="conda activate sch_man_nwinvasion-sweepfinder2; SweepFinder2 -lg 1000 $POP/$CHR"_"$POP.in $POP/$POP.sfs $POP/$CHR"_"$POP.sw2out"

        echo $CMD | qsub -V -cwd -S /bin/bash -q all.q -j y -N $CHR"_"$POP"_"sweep -o $POP/$CHR.log -pe smp 12
    done
done





