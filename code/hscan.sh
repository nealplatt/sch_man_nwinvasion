conda activate sch_man_nwinvasion-hscan

CONDA="conda activate sch_man_nwinvasion-hscan;"
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -pe smp 2 "    


cd ~/sch_man_nwinvasion

mkdir -p results/hscan/logs

#REAL DATA
for POP in new_world west_africa east_africa ; do
    mkdir results/hscan/$POP
        for i in $(seq 1 7 ); do
            CHROM="SM_V7_"$i
    
            #get chrom vcf
            vcftools \
                --vcf results/variant_filtration/smv7_ex_autosomes.vcf \
                --keep results/lists/$POP.list \
                --chr $CHROM \
                --recode \
                --recode-INFO-all \
                --stdout \
                >results/hscan/$POP/"$CHROM"_"$POP".vcf

            #convert vcf to hscan
            python code/vcf2hscan.py \
                results/hscan/$POP/"$CHROM"_"$POP".vcf \
                results/hscan/$POP/"$CHROM"_"$POP".hscan-in

            #convert to hscan
            bin/H-scan \
                -i results/hscan/$POP/"$CHROM"_"$POP".hscan-in \
                -g 10000 \
                >results/hscan/$POP/"$CHROM"_"$POP".hscan-out

    done
done


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     NOW RUN MSPRIME.SH TO GET NEUTRAL
#
#    these values are used to validate H
#    are scores comporable to actual data
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#now convert all vcfs to hscan fmt
for POP in new_world east_africa west_africa; do

    mkdir results/hscan/$POP"-msprime"

    for VCF in $(ls results/msprime/$POP/*probed.vcf); do
        BASE=$(basename $VCF)
        OUT=results/hscan/$POP"-msprime"/$(basename $VCF).hscan-in
        VCF2HSCAN_CMD="python code/vcf2hscan.py $VCF $OUT"
        
        echo "$CONDA $VCF2HSCAN_CMD" | $QSUB -N hscan2vcf_$BASE -o results/hscan/logs/"$BASE"_convert.log
    
    done
done

#batch submit hscan
for POP in new_world west_africa east_africa ; do

    for IN in $(ls results/hscan/$POP"-msprime"/*.hscan-in); do
        BASE=$(basename $IN .hscan-in)
        OUT=results/hscan/"$POP"-msprime/"$BASE".hscan-out
        
        bin/H-scan -i $IN -g 10000 >$OUT

    done
done


