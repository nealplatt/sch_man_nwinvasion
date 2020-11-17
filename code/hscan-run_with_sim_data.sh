conda activate sch_man_nwinvasion-hscan

CONDA="conda activate sch_man_nwinvasion-hscan"
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 1 "    


cd ~/sch_man_nwinvasion

mkdir -p results/hscan/logs

#convert all vcfs to hscan fmt
for POP in brazil tanzania senegal niger; do

    mkdir results/hscan/$POP"-sim"

    for VCF in $(ls results/msprime/$POP/*probed.vcf); do
        BASE=$(basename $VCF)
        OUT=results/hscan/$POP"-sim"/$(basename $VCF).hscan-in
        VCF2HSCAN_CMD="python code/vcf2hscan.py $VCF $OUT"
        
        echo "$CONDA; $VCF2HSCAN_CMD" | $QSUB -N hscan2vcf_$BASE -o results/hscan/logs/"$BASE"_convert.log
    
    done
done

#batch submit hscan
for POP in brazil tanzania senegal niger; do

    for IN in $(ls results/hscan/$POP"-sim"/*.hscan-in); do
        BASE=$(basename $IN .hscan-in)
        OUT=results/hscan/"$POP"-sim/"$BASE".hscan-out
        
        HSCAN_CMD="bin/H-scan -i $IN -g 10000 >$OUT"
        
        echo "$CONDA; $HSCAN_CMD" | $QSUB -N hscan_$BASE -o results/hscan/logs/"$BASE"_hscan.log

    done
done


